#include "bgzf.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

namespace {

enum class Fault {
    None,
    InflateInit,
    InflateEnd,
    DeflateInit,
    Deflate,
    DeflateEnd,
    OversizeLarge,
    OversizeAll
};

Fault fault = Fault::None;

class FailingBuffer : public std::streambuf {
protected:
    std::streamsize xsputn(const char *, std::streamsize) override {
        return 0;
    }

    int overflow(int) override {
        return traits_type::eof();
    }
};

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void expectFailure(const std::string &message,
                   const std::function<void()> &operation) {
    try {
        operation();
    } catch (const std::exception &error) {
        require(std::string(error.what()).find(message) != std::string::npos,
                "unexpected failure: " + std::string(error.what()));
        return;
    }
    throw std::runtime_error("expected failure: " + message);
}

std::string makeBgzf(const std::string &data) {
    std::ostringstream output(std::ios::binary);
    BgzfWriter writer(output);
    writer.write(reinterpret_cast<const uint8_t *>(data.data()), data.size());
    writer.finish();
    writer.finish();
    return output.str();
}

void readBgzf(const std::string &encoded, const std::string &expected) {
    std::istringstream input(encoded, std::ios::binary);
    BgzfReader reader(input);
    std::vector<uint8_t> decoded(expected.size());
    reader.readExact(decoded.data(), decoded.size(), "fixture");
    require(std::string(decoded.begin(), decoded.end()) == expected,
            "BGZF round trip failed");
    require(reader.read(decoded.data(), 1) == 0, "BGZF EOF was not reached");
    require(!reader.eofBlockMissing(), "BGZF EOF marker was not recognized");
}

void testZlibFailures(const std::string &encoded) {
    fault = Fault::InflateInit;
    expectFailure("initialize BGZF decompressor", [&]() {
        readBgzf(encoded, "teloscope");
    });

    fault = Fault::InflateEnd;
    expectFailure("finalize BGZF decompressor", [&]() {
        readBgzf(encoded, "teloscope");
    });

    fault = Fault::DeflateInit;
    expectFailure("initialize BGZF compressor", []() {
        makeBgzf("teloscope");
    });

    fault = Fault::Deflate;
    expectFailure("compress BGZF block", []() {
        makeBgzf("teloscope");
    });

    fault = Fault::DeflateEnd;
    expectFailure("finalize BGZF compressor", []() {
        makeBgzf("teloscope");
    });
}

void testOversizeFallback() {
    fault = Fault::OversizeLarge;
    const std::string encoded = makeBgzf("teloscope");
    fault = Fault::None;
    readBgzf(encoded, "teloscope");

    fault = Fault::OversizeAll;
    expectFailure("cannot fit compressed data", []() {
        makeBgzf("x");
    });
}

void testOutputFailures() {
    fault = Fault::None;
    FailingBuffer blockBuffer;
    std::ostream blockOutput(&blockBuffer);
    expectFailure("writing BAM output", [&]() {
        BgzfWriter writer(blockOutput);
        std::vector<uint8_t> data(65280, 0);
        writer.write(data.data(), data.size());
    });

    FailingBuffer finishBuffer;
    std::ostream finishOutput(&finishBuffer);
    expectFailure("finalizing BAM output", [&]() {
        BgzfWriter writer(finishOutput);
        writer.finish();
    });
}

} // namespace

extern "C" {

int __real_inflateInit2_(z_streamp, int, const char *, int);
int __real_inflateEnd(z_streamp);
int __real_deflateInit2_(z_streamp, int, int, int, int, int,
                         const char *, int);
int __real_deflate(z_streamp, int);
int __real_deflateEnd(z_streamp);
uLong __real_deflateBound(z_streamp, uLong);

int __wrap_inflateInit2_(z_streamp stream, int windowBits,
                         const char *version, int streamSize) {
    if (fault == Fault::InflateInit) return Z_MEM_ERROR;
    return __real_inflateInit2_(stream, windowBits, version, streamSize);
}

int __wrap_inflateEnd(z_streamp stream) {
    const int status = __real_inflateEnd(stream);
    return fault == Fault::InflateEnd ? Z_DATA_ERROR : status;
}

int __wrap_deflateInit2_(z_streamp stream, int level, int method,
                         int windowBits, int memLevel, int strategy,
                         const char *version, int streamSize) {
    if (fault == Fault::DeflateInit) return Z_MEM_ERROR;
    return __real_deflateInit2_(stream, level, method, windowBits, memLevel,
                                strategy, version, streamSize);
}

int __wrap_deflate(z_streamp stream, int flush) {
    if (fault == Fault::Deflate) return Z_STREAM_ERROR;
    const bool oversize = fault == Fault::OversizeAll ||
                          (fault == Fault::OversizeLarge && stream->avail_in > 1);
    if (!oversize) return __real_deflate(stream, flush);
    constexpr uInt outputSize = 66000;
    std::memset(stream->next_out, 0, outputSize);
    stream->next_in += stream->avail_in;
    stream->total_in += stream->avail_in;
    stream->avail_in = 0;
    stream->next_out += outputSize;
    stream->avail_out -= outputSize;
    stream->total_out += outputSize;
    return Z_STREAM_END;
}

int __wrap_deflateEnd(z_streamp stream) {
    const int status = __real_deflateEnd(stream);
    if (fault == Fault::OversizeLarge || fault == Fault::OversizeAll) {
        return Z_OK;
    }
    return fault == Fault::DeflateEnd ? Z_DATA_ERROR : status;
}

uLong __wrap_deflateBound(z_streamp stream, uLong size) {
    if (fault == Fault::OversizeLarge || fault == Fault::OversizeAll) {
        return 70000;
    }
    return __real_deflateBound(stream, size);
}

} // extern "C"

int main() {
    try {
        const std::string encoded = makeBgzf("teloscope");
        readBgzf(encoded, "teloscope");
        testZlibFailures(encoded);
        testOversizeFallback();
        testOutputFailures();
        std::cout << "PASS BGZF fault injection\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAIL BGZF fault injection: " << error.what() << '\n';
        return 1;
    }
}
