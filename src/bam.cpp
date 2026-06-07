#include "main.h"
#include "bam.h"

#include <algorithm>
#include <array>
#include <cinttypes>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "bgzf.h"
#include "global.h"
#include "read-filter.h"
#include "threadpool.h"

#ifdef _WIN32
#include <fcntl.h>
#endif

namespace {

constexpr size_t BAM_CORE_SIZE = 32;
constexpr size_t BAM_MAX_HEADER_TEXT = 1ULL << 30;
constexpr size_t BAM_MAX_REFERENCE_NAME = 1ULL << 20;
constexpr size_t BAM_MAX_RECORD_SIZE = 256ULL << 20;
constexpr size_t BAM_BATCH_BYTES = 32ULL << 20;

uint16_t getU16(const uint8_t *data) {
    return static_cast<uint16_t>(data[0]) |
           (static_cast<uint16_t>(data[1]) << 8);
}

uint32_t getU32(const uint8_t *data) {
    return static_cast<uint32_t>(data[0]) |
           (static_cast<uint32_t>(data[1]) << 8) |
           (static_cast<uint32_t>(data[2]) << 16) |
           (static_cast<uint32_t>(data[3]) << 24);
}

int32_t getI32(const uint8_t *data) {
    uint32_t value = getU32(data);
    int32_t result;
    std::memcpy(&result, &value, sizeof(result));
    return result;
}

struct BamRecord {
    std::vector<uint8_t> raw;
    std::string sequence;
    bool hasSequence = false;
};

struct BamSubsetStats {
    uint64_t totalRecords = 0;
    uint64_t passedRecords = 0;
    uint64_t missingSequenceRecords = 0;
    bool missingEofBlock = false;
};

void copyBytes(BgzfReader &reader, BgzfWriter &writer, size_t size,
               const char *context) {
    std::array<uint8_t, 65536> buffer{};
    while (size > 0) {
        const size_t chunk = std::min(size, buffer.size());
        reader.readExact(buffer.data(), chunk, context);
        writer.write(buffer.data(), chunk);
        size -= chunk;
    }
}

void copyBamHeader(BgzfReader &reader, BgzfWriter &writer) {
    std::array<uint8_t, 4> word{};
    reader.readExact(word.data(), word.size(), "magic");
    if (word != std::array<uint8_t, 4>{'B', 'A', 'M', 1}) {
        throw std::runtime_error("invalid BAM magic");
    }
    writer.write(word.data(), word.size());

    reader.readExact(word.data(), word.size(), "header length");
    const int32_t textLength = getI32(word.data());
    if (textLength < 0 || static_cast<size_t>(textLength) > BAM_MAX_HEADER_TEXT) {
        throw std::runtime_error("invalid BAM header text length");
    }
    writer.write(word.data(), word.size());
    copyBytes(reader, writer, static_cast<size_t>(textLength), "header text");

    reader.readExact(word.data(), word.size(), "reference count");
    const int32_t referenceCount = getI32(word.data());
    if (referenceCount < 0) {
        throw std::runtime_error("invalid BAM reference count");
    }
    writer.write(word.data(), word.size());

    for (int32_t i = 0; i < referenceCount; ++i) {
        reader.readExact(word.data(), word.size(), "reference name length");
        const int32_t nameLength = getI32(word.data());
        if (nameLength <= 0 ||
            static_cast<size_t>(nameLength) > BAM_MAX_REFERENCE_NAME) {
            throw std::runtime_error("invalid BAM reference name length");
        }
        writer.write(word.data(), word.size());

        std::vector<uint8_t> name(static_cast<size_t>(nameLength));
        reader.readExact(name.data(), name.size(), "reference name");
        if (name.back() != 0) {
            throw std::runtime_error("BAM reference name is not NUL-terminated");
        }
        writer.write(name.data(), name.size());

        reader.readExact(word.data(), word.size(), "reference length");
        if (getI32(word.data()) < 0) {
            throw std::runtime_error("invalid BAM reference length");
        }
        writer.write(word.data(), word.size());
    }
}

std::string decodeSequence(const std::vector<uint8_t> &raw) {
    const uint8_t *core = raw.data() + 4;
    const size_t payloadSize = raw.size() - 4;
    const size_t readNameLength = core[8];
    const size_t cigarCount = getU16(core + 12);
    const int32_t signedSequenceLength = getI32(core + 16);
    if (readNameLength == 0 || signedSequenceLength < 0) {
        throw std::runtime_error("invalid BAM record lengths");
    }
    const size_t sequenceLength = static_cast<size_t>(signedSequenceLength);

    size_t sequenceOffset = BAM_CORE_SIZE;
    if (readNameLength > payloadSize - sequenceOffset) {
        throw std::runtime_error("BAM read name exceeds block_size");
    }
    sequenceOffset += readNameLength;
    const size_t cigarBytes = cigarCount * 4;
    if (cigarBytes > payloadSize - sequenceOffset) {
        throw std::runtime_error("BAM CIGAR exceeds block_size");
    }
    sequenceOffset += cigarBytes;
    const size_t packedLength = sequenceLength / 2 + sequenceLength % 2;
    if (packedLength > payloadSize - sequenceOffset) {
        throw std::runtime_error("BAM sequence exceeds block_size");
    }
    const size_t qualityOffset = sequenceOffset + packedLength;
    if (sequenceLength > payloadSize - qualityOffset) {
        throw std::runtime_error("BAM record fields exceed block_size");
    }
    if (core[BAM_CORE_SIZE + readNameLength - 1] != 0) {
        throw std::runtime_error("BAM read name is not NUL-terminated");
    }

    static constexpr char bases[] = "=ACMGRSVTWYHKDBN";
    std::string sequence(sequenceLength, 'N');
    for (size_t i = 0; i < sequenceLength; ++i) {
        const uint8_t packed = core[sequenceOffset + i / 2];
        const uint8_t code = (i % 2 == 0) ? (packed >> 4) : (packed & 0x0f);
        sequence[i] = bases[code];
    }
    return sequence;
}

bool readRecord(BgzfReader &reader, BamRecord &record) {
    std::array<uint8_t, 4> sizeBytes{};
    const size_t got = reader.read(sizeBytes.data(), sizeBytes.size());
    if (got == 0) return false;
    if (got != sizeBytes.size()) {
        throw std::runtime_error("truncated BAM record size");
    }

    const int32_t blockSize = getI32(sizeBytes.data());
    if (blockSize < static_cast<int32_t>(BAM_CORE_SIZE) ||
        static_cast<size_t>(blockSize) > BAM_MAX_RECORD_SIZE) {
        throw std::runtime_error("invalid BAM record block_size");
    }

    record.raw.resize(static_cast<size_t>(blockSize) + sizeBytes.size());
    std::copy(sizeBytes.begin(), sizeBytes.end(), record.raw.begin());
    reader.readExact(record.raw.data() + sizeBytes.size(),
                     static_cast<size_t>(blockSize), "record");
    record.sequence = decodeSequence(record.raw);
    record.hasSequence = !record.sequence.empty();
    return true;
}

BamSubsetStats subsetBam(std::istream &input,
                         std::ostream &output,
                         const UserInputTeloscope &userInput) {
    BgzfReader reader(input);
    BgzfWriter writer(output);
    copyBamHeader(reader, writer);

    const uint32_t threads = std::max<uint32_t>(1, threadPool.totalThreads());
    const size_t recordsPerBatch = std::min<size_t>(
        2048, std::max<size_t>(256, static_cast<size_t>(threads) * 32));

    std::vector<BamRecord> batch;
    batch.reserve(recordsPerBatch);
    size_t batchBytes = 0;
    BamSubsetStats stats;

    auto processBatch = [&]() {
        if (batch.empty()) return;

        std::vector<uint8_t> passed(batch.size(), 0);
        const size_t chunkCount = std::min<size_t>(threads, batch.size());
        const size_t chunkSize = (batch.size() + chunkCount - 1) / chunkCount;

        for (size_t chunk = 0; chunk < chunkCount; ++chunk) {
            const size_t start = chunk * chunkSize;
            const size_t end = std::min(batch.size(), start + chunkSize);
            if (start >= end) continue;

            threadPool.queueJob([&, start, end]() {
                ReadTelomereFilter filter(userInput);
                for (size_t i = start; i < end; ++i) {
                    if (batch[i].hasSequence &&
                        filter.matches(std::move(batch[i].sequence))) {
                        passed[i] = 1;
                    }
                }
                return true;
            });
        }
        jobWait(threadPool);

        for (size_t i = 0; i < batch.size(); ++i) {
            stats.totalRecords++;
            if (!batch[i].hasSequence) {
                stats.missingSequenceRecords++;
            } else if (passed[i]) {
                writer.write(batch[i].raw.data(), batch[i].raw.size());
                stats.passedRecords++;
            }
        }

        batch.clear();
        batchBytes = 0;
    };

    BamRecord record;
    while (readRecord(reader, record)) {
        if (!batch.empty() &&
            (batch.size() >= recordsPerBatch ||
             record.raw.size() > BAM_BATCH_BYTES - std::min(batchBytes, BAM_BATCH_BYTES))) {
            processBatch();
        }
        batchBytes += record.raw.size();
        batch.push_back(std::move(record));
        record = BamRecord();
    }
    processBatch();
    writer.finish();

    stats.missingEofBlock = reader.eofBlockMissing();
    return stats;
}

} // namespace

void runBamSubsetMode(const UserInputTeloscope &userInput) {
#ifdef _WIN32
    _setmode(_fileno(stdin), _O_BINARY);
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::ifstream inputFile;
    std::istream *input = &std::cin;
    if (!userInput.inSequence.empty()) {
        inputFile.open(userInput.inSequence, std::ios::binary);
        if (!inputFile.is_open()) {
            throw std::runtime_error("cannot open BAM input '" + userInput.inSequence + "'");
        }
        input = &inputFile;
    }

    std::ofstream outputFile;
    std::ostream *output = &std::cout;
    std::string outputPath;
    if (userInput.outRouteSet) {
        const std::filesystem::path inputName(userInput.inSequenceName);
        outputPath = userInput.outRoute + "/" + inputName.stem().string() + "_telomeric.bam";
        outputFile.open(outputPath, std::ios::binary);
        if (!outputFile.is_open()) {
            throw std::runtime_error("cannot write telomeric records to '" + outputPath + "'");
        }
        output = &outputFile;
    }

    BamSubsetStats stats;
    try {
        stats = subsetBam(*input, *output, userInput);
    } catch (...) {
        if (outputFile.is_open()) {
            outputFile.close();
            std::error_code removeError;
            std::filesystem::remove(outputPath, removeError);
        }
        throw;
    }

    if (stats.missingEofBlock) {
        fprintf(stderr, "Warning: BAM input is missing the BGZF EOF marker.\n");
    }
    if (stats.missingSequenceRecords > 0) {
        fprintf(stderr, "BAM subset: skipped %" PRIu64 " record%s without SEQ.\n",
                stats.missingSequenceRecords,
                stats.missingSequenceRecords == 1 ? "" : "s");
    }
    fprintf(stderr, "BAM subset: kept %" PRIu64 " of %" PRIu64 " records.\n",
            stats.passedRecords, stats.totalRecords);
    if (!outputPath.empty()) {
        fprintf(stderr, "Wrote telomeric records to %s.\n", outputPath.c_str());
    }
}
