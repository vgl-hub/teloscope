#include "bgzf.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <exception>
#include <istream>
#include <mutex>
#include <ostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

#include <zlib.h>

namespace {

constexpr size_t MAX_BLOCK_SIZE = 65536;
constexpr size_t MAX_UNCOMPRESSED = 65280;
constexpr size_t GROUP_BLOCKS = 64;
constexpr std::array<uint8_t, 28> EOF_BLOCK = {
    31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67,
    2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

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

void putU16(uint8_t *data, uint16_t value) {
    data[0] = static_cast<uint8_t>(value);
    data[1] = static_cast<uint8_t>(value >> 8);
}

void putU32(uint8_t *data, uint32_t value) {
    data[0] = static_cast<uint8_t>(value);
    data[1] = static_cast<uint8_t>(value >> 8);
    data[2] = static_cast<uint8_t>(value >> 16);
    data[3] = static_cast<uint8_t>(value >> 24);
}

// raw blocks filled by the reader thread, inflated by one worker, drained in file order
struct Group {
    bool ready = false;
    std::vector<uint8_t> raw;
    std::vector<size_t> ends;
    std::vector<uint8_t> out;
    std::exception_ptr error;
};

void inflateBlock(z_stream &stream, const uint8_t *block, size_t totalSize,
                  std::vector<uint8_t> &out) {
    const uint8_t flags = block[3];
    size_t payloadOffset = 12 + getU16(block + 10);
    const size_t footerOffset = totalSize - 8;
    auto skipZeroTerminated = [&](const char *fieldName) {
        while (payloadOffset < footerOffset && block[payloadOffset] != 0) {
            payloadOffset++;
        }
        if (payloadOffset == footerOffset) {
            throw std::runtime_error(std::string("unterminated BGZF ") + fieldName);
        }
        payloadOffset++;
    };
    if (flags & 0x08) skipZeroTerminated("filename");
    if (flags & 0x10) skipZeroTerminated("comment");
    if (flags & 0x02) {
        if (footerOffset - payloadOffset < 2) {
            throw std::runtime_error("truncated BGZF header checksum");
        }
        const uint16_t expectedHeaderCrc = getU16(block + payloadOffset);
        const uLong actualHeaderCrc = crc32(
            crc32(0L, Z_NULL, 0), block, static_cast<uInt>(payloadOffset));
        if (static_cast<uint16_t>(actualHeaderCrc) != expectedHeaderCrc) {
            throw std::runtime_error("BGZF header checksum mismatch");
        }
        payloadOffset += 2;
    }
    const uint32_t expectedCrc = getU32(block + footerOffset);
    const uint32_t expectedSize = getU32(block + footerOffset + 4);
    if (expectedSize > MAX_BLOCK_SIZE) {
        throw std::runtime_error("BGZF uncompressed block is too large");
    }

    const size_t start = out.size();
    out.resize(start + std::max<uint32_t>(expectedSize, 1));
    inflateReset(&stream);
    stream.next_in = const_cast<Bytef *>(block + payloadOffset);
    stream.avail_in = static_cast<uInt>(footerOffset - payloadOffset);
    stream.next_out = out.data() + start;
    stream.avail_out = static_cast<uInt>(out.size() - start);
    if (inflate(&stream, Z_FINISH) != Z_STREAM_END ||
        stream.total_out != expectedSize ||
        stream.total_in != footerOffset - payloadOffset) {
        throw std::runtime_error("invalid BGZF deflate payload");
    }
    const uLong actualCrc = crc32(
        crc32(0L, Z_NULL, 0),
        expectedSize == 0 ? Z_NULL : out.data() + start,
        expectedSize);
    if (actualCrc != expectedCrc) {
        throw std::runtime_error("BGZF checksum mismatch");
    }
    out.resize(start + expectedSize);
}

void inflateGroup(Group &group) {
    group.out.clear();
    z_stream stream{};
    if (inflateInit2(&stream, -15) != Z_OK) {
        throw std::runtime_error("could not initialize BGZF decompressor");
    }
    size_t start = 0;
    size_t outStart = 0;
    try {
        for (const size_t end : group.ends) {
            outStart = group.out.size();
            inflateBlock(stream, group.raw.data() + start, end - start, group.out);
            start = end;
        }
    } catch (...) {
        inflateEnd(&stream);
        group.out.resize(outStart);
        throw;
    }
    if (inflateEnd(&stream) != Z_OK) {
        throw std::runtime_error("could not finalize BGZF decompressor");
    }
}

} // namespace

// own threads rather than the global pool, whose jobWait the consumer already blocks on
struct BgzfReader::Pipeline {
    std::istream &input;
    std::vector<Group> groups;
    std::mutex mutex;
    std::condition_variable freed, filled, ready;
    std::deque<size_t> queue;
    size_t produced = 0;
    size_t consumed = 0;
    std::exception_ptr pendingError;
    bool done = false;
    std::exception_ptr readerError;
    bool sawEofBlock = false;
    std::atomic<bool> stop{false};
    std::thread reader;
    std::vector<std::thread> workers;

    Pipeline(std::istream &stream, unsigned threads);
    ~Pipeline();
    void readPhysical(uint8_t *data, size_t size, const char *context);
    bool readBlock(Group &group);
    void produce();
    void work();
};

BgzfReader::Pipeline::Pipeline(std::istream &stream, unsigned threads)
    : input(stream), groups(2 * std::max(1u, threads)) {
    reader = std::thread(&Pipeline::produce, this);
    for (unsigned i = 0; i < std::max(1u, threads); ++i) {
        workers.emplace_back(&Pipeline::work, this);
    }
}

BgzfReader::Pipeline::~Pipeline() {
    {
        std::lock_guard<std::mutex> lock(mutex);
        stop = true;
    }
    freed.notify_all();
    filled.notify_all();
    reader.join();
    for (std::thread &worker : workers) worker.join();
}

void BgzfReader::Pipeline::readPhysical(uint8_t *data, size_t size, const char *context) {
    input.read(reinterpret_cast<char *>(data), static_cast<std::streamsize>(size));
    if (static_cast<size_t>(input.gcount()) != size) {
        throw std::runtime_error(std::string("truncated BGZF ") + context);
    }
}

bool BgzfReader::Pipeline::readBlock(Group &group) {
    if (stop) return false;
    std::array<uint8_t, 12> prefix{};
    input.read(reinterpret_cast<char *>(prefix.data()),
               static_cast<std::streamsize>(prefix.size()));
    const size_t got = static_cast<size_t>(input.gcount());
    if (got == 0 && input.eof()) return false;
    if (got != prefix.size()) {
        throw std::runtime_error("truncated BGZF header");
    }
    if (prefix[0] != 31 || prefix[1] != 139 || prefix[2] != 8) {
        throw std::runtime_error("input is not BGZF-compressed BAM");
    }

    const uint8_t flags = prefix[3];
    if ((flags & 0x04) == 0 || (flags & 0xe0) != 0) {
        throw std::runtime_error("invalid BGZF gzip flags");
    }

    const uint16_t extraLength = getU16(prefix.data() + 10);
    const size_t start = group.raw.size();
    const size_t bytesRead = prefix.size() + extraLength;
    group.raw.resize(start + bytesRead);
    std::copy(prefix.begin(), prefix.end(), group.raw.begin() + start);
    uint8_t *extra = group.raw.data() + start + prefix.size();
    if (extraLength > 0) {
        readPhysical(extra, extraLength, "extra field");
    }

    bool foundBlockSize = false;
    uint16_t storedBlockSize = 0;
    for (size_t pos = 0; pos < extraLength;) {
        if (extraLength - pos < 4) {
            throw std::runtime_error("malformed BGZF extra field");
        }
        const uint16_t subfieldLength = getU16(extra + pos + 2);
        const size_t end = pos + 4 + subfieldLength;
        if (end > extraLength) {
            throw std::runtime_error("malformed BGZF extra subfield");
        }
        if (extra[pos] == 'B' && extra[pos + 1] == 'C') {
            if (subfieldLength != 2 || foundBlockSize) {
                throw std::runtime_error("invalid BGZF BC subfield");
            }
            storedBlockSize = getU16(extra + pos + 4);
            foundBlockSize = true;
        }
        pos = end;
    }
    if (!foundBlockSize) {
        throw std::runtime_error("BGZF block is missing the BC subfield");
    }

    const size_t totalSize = static_cast<size_t>(storedBlockSize) + 1;
    if (totalSize > MAX_BLOCK_SIZE || totalSize < bytesRead + 8) {
        throw std::runtime_error("invalid BGZF block size");
    }
    group.raw.resize(start + totalSize);
    readPhysical(group.raw.data() + start + bytesRead, totalSize - bytesRead, "block");
    group.ends.push_back(group.raw.size());
    if (totalSize == EOF_BLOCK.size() &&
        std::equal(EOF_BLOCK.begin(), EOF_BLOCK.end(), group.raw.begin() + start)) {
        sawEofBlock = true;
    }
    return true;
}

void BgzfReader::Pipeline::produce() {
    bool eof = false;
    while (!eof) {
        Group &group = groups[produced % groups.size()];
        std::unique_lock<std::mutex> lock(mutex);
        freed.wait(lock, [&] { return stop || produced < consumed + groups.size(); });
        if (stop) return;
        lock.unlock();
        group.raw.clear();
        group.ends.clear();
        try {
            while (!eof && group.ends.size() < GROUP_BLOCKS) eof = !readBlock(group);
        } catch (...) {
            readerError = std::current_exception();
            eof = true;
        }
        lock.lock();
        if (!group.ends.empty()) {
            queue.push_back(produced++);
            filled.notify_one();
        }
        if (eof) {
            done = true;
            ready.notify_one();
        }
    }
}

void BgzfReader::Pipeline::work() {
    while (true) {
        std::unique_lock<std::mutex> lock(mutex);
        filled.wait(lock, [&] { return stop || !queue.empty(); });
        if (stop) return;
        Group &group = groups[queue.front() % groups.size()];
        queue.pop_front();
        lock.unlock();
        try {
            inflateGroup(group);
        } catch (...) {
            group.error = std::current_exception();
        }
        lock.lock();
        group.ready = true;
        ready.notify_one();
    }
}

BgzfReader::BgzfReader(std::istream &stream, unsigned threads)
    : pipeline(new Pipeline(stream, threads)) {}

BgzfReader::~BgzfReader() = default;

bool BgzfReader::loadGroup() {
    Pipeline &p = *pipeline;
    if (p.pendingError) std::rethrow_exception(p.pendingError);
    Group &group = p.groups[p.consumed % p.groups.size()];
    std::unique_lock<std::mutex> lock(p.mutex);
    p.ready.wait(lock, [&] { return group.ready || (p.done && p.consumed == p.produced); });
    if (!group.ready) {
        if (p.readerError) std::rethrow_exception(p.readerError);
        physicalEof = true;
        sawEofBlock = p.sawEofBlock;
        return false;
    }
    // a bad block's error follows the good blocks before it, as the serial reader did
    p.pendingError = group.error;
    std::swap(uncompressed, group.out);
    offset = 0;
    group.ready = false;
    p.consumed++;
    p.freed.notify_one();
    return true;
}

size_t BgzfReader::read(uint8_t *data, size_t size) {
    size_t copied = 0;
    while (copied < size) {
        if (offset == uncompressed.size() && !loadGroup()) {
            break;
        }
        const size_t available = uncompressed.size() - offset;
        const size_t take = std::min(size - copied, available);
        std::copy_n(uncompressed.data() + offset, take, data + copied);
        offset += take;
        copied += take;
    }
    return copied;
}

void BgzfReader::readExact(uint8_t *data, size_t size, const char *context) {
    if (read(data, size) != size) {
        throw std::runtime_error(std::string("truncated BAM ") + context);
    }
}

bool BgzfReader::eofBlockMissing() const {
    return physicalEof && !sawEofBlock;
}

BgzfWriter::BgzfWriter(std::ostream &stream) : output(stream) {
    pending.reserve(MAX_UNCOMPRESSED);
}

std::vector<uint8_t> BgzfWriter::compressRaw(const uint8_t *data, size_t size) {
    z_stream stream{};
    if (deflateInit2(&stream, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
                     -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        throw std::runtime_error("could not initialize BGZF compressor");
    }

    std::vector<uint8_t> compressed(deflateBound(&stream, size));
    stream.next_in = const_cast<Bytef *>(data);
    stream.avail_in = static_cast<uInt>(size);
    stream.next_out = compressed.data();
    stream.avail_out = static_cast<uInt>(compressed.size());
    const int status = deflate(&stream, Z_FINISH);
    const size_t compressedSize = stream.total_out;
    const int endStatus = deflateEnd(&stream);
    if (status != Z_STREAM_END) {
        throw std::runtime_error("could not compress BGZF block");
    }
    if (endStatus != Z_OK) {
        throw std::runtime_error("could not finalize BGZF compressor");
    }
    compressed.resize(compressedSize);
    return compressed;
}

void BgzfWriter::emitBlock(const uint8_t *data, size_t size) {
    std::vector<uint8_t> compressed = compressRaw(data, size);
    const size_t totalSize = 18 + compressed.size() + 8;
    if (totalSize > MAX_BLOCK_SIZE) {
        if (size <= 1) {
            throw std::runtime_error("BGZF block cannot fit compressed data");
        }
        const size_t split = size / 2;
        emitBlock(data, split);
        emitBlock(data + split, size - split);
        return;
    }

    std::vector<uint8_t> block(totalSize, 0);
    block[0] = 31;
    block[1] = 139;
    block[2] = 8;
    block[3] = 4;
    block[9] = 255;
    putU16(block.data() + 10, 6);
    block[12] = 'B';
    block[13] = 'C';
    putU16(block.data() + 14, 2);
    putU16(block.data() + 16, static_cast<uint16_t>(totalSize - 1));
    std::copy(compressed.begin(), compressed.end(), block.begin() + 18);

    const uLong crc = crc32(crc32(0L, Z_NULL, 0), data, static_cast<uInt>(size));
    putU32(block.data() + totalSize - 8, static_cast<uint32_t>(crc));
    putU32(block.data() + totalSize - 4, static_cast<uint32_t>(size));
    output.write(reinterpret_cast<const char *>(block.data()),
                 static_cast<std::streamsize>(block.size()));
    if (!output.good()) {
        throw std::runtime_error("failed while writing BAM output");
    }
}

void BgzfWriter::flushPending() {
    if (!pending.empty()) {
        emitBlock(pending.data(), pending.size());
        pending.clear();
    }
}

void BgzfWriter::write(const uint8_t *data, size_t size) {
    while (size > 0) {
        const size_t capacity = MAX_UNCOMPRESSED - pending.size();
        const size_t take = std::min(size, capacity);
        pending.insert(pending.end(), data, data + take);
        data += take;
        size -= take;
        if (pending.size() == MAX_UNCOMPRESSED) {
            flushPending();
        }
    }
}

void BgzfWriter::finish() {
    if (finished) return;
    flushPending();
    output.write(reinterpret_cast<const char *>(EOF_BLOCK.data()),
                 static_cast<std::streamsize>(EOF_BLOCK.size()));
    output.flush();
    if (!output.good()) {
        throw std::runtime_error("failed while finalizing BAM output");
    }
    finished = true;
}
