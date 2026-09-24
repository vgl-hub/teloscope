#ifndef BGZF_H
#define BGZF_H

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <memory>
#include <vector>

class BgzfReader {
    struct Pipeline;
    std::unique_ptr<Pipeline> pipeline;
    std::vector<uint8_t> uncompressed;
    size_t offset = 0;
    bool physicalEof = false;
    bool sawEofBlock = false;

    bool loadGroup();

public:
    // threads counts inflate workers; bytes are still delivered in file order
    explicit BgzfReader(std::istream &stream, unsigned threads = 1);
    ~BgzfReader();
    size_t read(uint8_t *data, size_t size);
    void readExact(uint8_t *data, size_t size, const char *context);
    bool eofBlockMissing() const;
};

class BgzfWriter {
    std::ostream &output;
    std::vector<uint8_t> pending;
    bool finished = false;

    std::vector<uint8_t> compressRaw(const uint8_t *data, size_t size);
    void emitBlock(const uint8_t *data, size_t size);
    void flushPending();

public:
    explicit BgzfWriter(std::ostream &stream);
    void write(const uint8_t *data, size_t size);
    void finish();
};

#endif /* BGZF_H */
