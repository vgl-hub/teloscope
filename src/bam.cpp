#include "main.h"
#include "bam.h"

#include <algorithm>
#include <array>
#include <atomic>
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
#include "teloscope.h"
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

// validates the record layout; IUPAC codes and '=' become N since the engine matches ACGT only
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

    static constexpr char bases[] = "NACNGNNNTNNNNNNN";
    std::string sequence(sequenceLength, 'N');
    for (size_t i = 0; i < sequenceLength; ++i) {
        const uint8_t packed = core[sequenceOffset + i / 2];
        const uint8_t code = (i % 2 == 0) ? (packed >> 4) : (packed & 0x0f);
        sequence[i] = bases[code];
    }
    return sequence;
}

// size and bounds are validated here, on the reader thread; SEQ is decoded in the worker jobs
bool readRawRecord(BgzfReader &reader, std::vector<uint8_t> &raw) {
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

    raw.resize(static_cast<size_t>(blockSize) + sizeBytes.size());
    std::copy(sizeBytes.begin(), sizeBytes.end(), raw.begin());
    reader.readExact(raw.data() + sizeBytes.size(), static_cast<size_t>(blockSize), "record");
    return true;
}

uint16_t recordFlag(const std::vector<uint8_t> &raw) {
    return getU16(raw.data() + 4 + 14);
}

std::string recordName(const std::vector<uint8_t> &raw) {
    const uint8_t *core = raw.data() + 4;
    return std::string(reinterpret_cast<const char *>(core + BAM_CORE_SIZE), core[8] - 1);
}

// only valid after decodeSequence accepted the record layout
bool hardClipped(const std::vector<uint8_t> &raw) {
    const uint8_t *core = raw.data() + 4;
    const size_t cigarCount = getU16(core + 12);
    if (cigarCount == 0) return false;
    const uint8_t *cigar = core + BAM_CORE_SIZE + core[8];
    return (getU32(cigar) & 0xf) == 5 || (getU32(cigar + 4 * (cigarCount - 1)) & 0xf) == 5;
}

// per record: measured (primary, SEQ present, no hard clip) or not; kept when a block exists or the keep scan passes
struct BamRecordResult {
    std::vector<TelomereBlock> blocks;
    bool hasSequence = false;
    bool measured = false;
    bool kept = false;
    bool malformed = false;
};

} // namespace

void readBamReads(const UserInputTeloscope &userInput, std::ostream &subset,
                  std::ostream &bed, ReadTlStats &stats) {
#ifdef _WIN32
    _setmode(_fileno(stdin), _O_BINARY);
#endif
    std::ifstream inputFile;
    std::istream *input = &std::cin;
    if (!userInput.inSequence.empty()) {
        inputFile.open(userInput.inSequence, std::ios::binary);
        input = &inputFile;
    }

    const uint32_t threads = std::max<uint32_t>(1, threadPool.totalThreads());
    BgzfReader reader(*input, threads);
    BgzfWriter writer(subset);
    copyBamHeader(reader, writer);

    const size_t recordsPerBatch = std::min<size_t>(
        2048, std::max<size_t>(256, static_cast<size_t>(threads) * 32));

    // two batches alternate: the reader fills one while the workers scan the other
    struct Batch {
        std::vector<std::vector<uint8_t>> records;
        std::vector<BamRecordResult> results;
        std::atomic<size_t> next{0};
    };
    Batch batches[2];
    std::vector<uint8_t> raw; // a record read ahead that did not fit the batch being filled
    uint64_t totalRecords = 0, missingSequence = 0, secondary = 0, clipped = 0;

    auto readBatch = [&](Batch &batch) {
        batch.records.clear();
        size_t bytes = 0;
        while (!raw.empty() || readRawRecord(reader, raw)) {
            if (!batch.records.empty() &&
                (batch.records.size() >= recordsPerBatch ||
                 raw.size() > BAM_BATCH_BYTES - std::min(bytes, BAM_BATCH_BYTES))) {
                return;
            }
            bytes += raw.size();
            batch.records.push_back(std::move(raw));
            raw.clear();
        }
    };

    // one job per thread; each pulls the next record so a run of long reads does not idle the rest
    auto scanBatch = [&](Batch &batch) {
        const size_t count = batch.records.size();
        batch.results.assign(count, {});
        batch.next = 0;
        for (size_t job = 0; job < std::min<size_t>(threads, count); ++job) {
            threadPool.queueJob([&, count]() {
                ReadTelomereScanner scanner(userInput);
                ReadTelomereFilter filter(userInput);
                for (size_t i = batch.next++; i < count; i = batch.next++) {
                    const std::vector<uint8_t> &record = batch.records[i];
                    BamRecordResult &result = batch.results[i];
                    std::string sequence;
                    try {
                        sequence = decodeSequence(record);
                    } catch (const std::exception &) {
                        result.malformed = true;
                        continue;
                    }
                    if (sequence.empty()) continue;
                    result.hasSequence = true;

                    const uint16_t flag = recordFlag(record);
                    result.measured = !(flag & 0x900) && !hardClipped(record);
                    if (result.measured) // reverse-strand records are measured in sequencing orientation
                        result.blocks = scanner.scan((flag & 0x10) ? revCom(sequence) : sequence);
                    result.kept = !result.blocks.empty() || filter.matches(std::move(sequence));
                }
                return true;
            });
        }
    };

    // written sequentially here, in input record order, regardless of -j
    auto writeBatch = [&](const Batch &batch) {
        for (size_t i = 0; i < batch.records.size(); ++i) {
            const std::vector<uint8_t> &record = batch.records[i];
            const BamRecordResult &result = batch.results[i];
            if (result.malformed) decodeSequence(record); // rethrows the record's error
            totalRecords++;
            if (!result.hasSequence) {
                missingSequence++;
                continue;
            }
            if (result.kept) {
                writer.write(record.data(), record.size());
                stats.readsKept++;
            }
            if (!result.measured) {
                if (recordFlag(record) & 0x900) secondary++;
                else clipped++;
                continue;
            }
            stats.readsMeasured++;
            const std::string name = recordName(record);
            const uint64_t readLen = static_cast<uint64_t>(getI32(record.data() + 4 + 16));
            for (const TelomereBlock &block : result.blocks) {
                writeReadTelomereRow(bed, userInput, name, readLen, block, stats);
            }
        }
    };

    // batch N+1 is read and batch N-1 written while the workers scan batch N
    Batch *scanning = nullptr;
    try {
        for (size_t fill = 0; ; fill ^= 1) {
            Batch &next = batches[fill];
            readBatch(next);
            if (scanning != nullptr) jobWait(threadPool);
            if (!next.records.empty()) scanBatch(next);
            if (scanning != nullptr) writeBatch(*scanning);
            if (next.records.empty()) break;
            scanning = &next;
        }
    } catch (...) {
        jobWait(threadPool); // the running jobs still use a batch on this frame
        throw;
    }
    writer.finish();

    if (reader.eofBlockMissing()) {
        fprintf(stderr, "Warning: BAM input is missing the BGZF EOF marker.\n");
    }
    if (missingSequence > 0) {
        fprintf(stderr, "BAM: skipped %" PRIu64 " record%s without SEQ.\n",
                missingSequence, missingSequence == 1 ? "" : "s");
    }
    if (secondary > 0) {
        fprintf(stderr, "BAM: %" PRIu64 " secondary/supplementary record%s not measured.\n",
                secondary, secondary == 1 ? "" : "s");
    }
    if (clipped > 0) {
        fprintf(stderr, "BAM: %" PRIu64 " hard-clipped record%s not measured.\n",
                clipped, clipped == 1 ? "" : "s");
    }
    fprintf(stderr, "BAM: kept %" PRIu64 " of %" PRIu64 " records.\n", stats.readsKept, totalRecords);
}
