#ifndef INPUT_H
#define INPUT_H

#include <algorithm>
#include <stdint.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>

#include "global.h"
#include "log.h"
#include "struct.h"

class InSequences;
struct ReadTlStats;

// set by content detection in main: FASTQ or BAM input is subset and measured, anything else is an assembly
enum class ReadInput : uint8_t { none, fastq, bam };

struct UserInputTeloscope : UserInput {
    
    std::string inSequencePrefix;
    std::string inSequenceName;
    std::string outRoute;
    bool outRouteSet = false;
    std::vector<std::string> includeBedFiles;
    std::vector<std::string> excludeBedFiles;
    std::vector<std::string> includePrefixes;
    std::vector<std::string> excludePrefixes;
    bool sequenceFilterActive = false;
    bool chrOnly = false;
    uint64_t filterInputCount = 0;
    uint64_t filterSelectedCount = 0;
    std::unordered_map<std::string, uint8_t> hammingDistances;

    std::string canonicalFwd = "CCCTAA";
    std::string canonicalRev = "TTAGGG";
    unsigned short int canonicalSize = 6;
    std::vector<std::string> rawPatterns = {"TTAGGG", "CCCTAA"};
    std::vector<std::string> patterns = {"TTAGGG", "CCCTAA"};
    std::vector<std::pair<std::string, bool>> patternInfo; // (pattern, isForward)

    uint32_t windowSize = 1000;
    uint32_t step = 1000;
    uint32_t terminalLimit = 50000;
    bool terminalLimitSet = false;
    uint8_t editDistance = 1;

    uint32_t maxMatchDist = 50;
    uint32_t minBlockLen = 300;
    uint32_t maxBlockDist = 1000;
    uint32_t minBlockCounts = 2;
    uint32_t terminalTolerance = 3000;
    bool terminalToleranceSet = false;
    float minBlockDensity = 0.5f;
    float labelThreshold = 0.667f;

    bool outFasta = false;
    bool outWinRepeats = false;
    bool outGC = false;
    bool outEntropy = false;
    bool outMatches = false;
    bool ultraFastMode = true;
    bool manualCuration = false;
    bool outPlotReport = false;
    ReadInput readInput = ReadInput::none;

    double maxMem = 0;
    std::string prefix = ".", outFile = "";
};

bool isGfaAssemblyPath(const std::string &path);

constexpr size_t READ_BATCH_BYTES = 32ULL << 20; // batch byte cap shared by the BAM and FASTQ readers

// enough records to occupy every worker thread while keeping batches small
inline size_t defaultRecordsPerBatch(uint32_t threads) {
    return std::min<size_t>(2048, std::max<size_t>(256, static_cast<size_t>(threads) * 32));
}

// batch N+1 is read and batch N-1 written while the workers scan batch N
template <typename Batch, typename Read, typename Scan, typename Write>
void scanBatchesDoubleBuffered(Batch (&batches)[2], Read readBatch, Scan scanBatch, Write writeBatch) {
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
}

class Input {

    UserInputTeloscope userInput;
    std::shared_ptr<std::istream> stream;

public:

    std::vector<Log> logs;

    void load(UserInputTeloscope userInput);
    void read(InSequences &inSequence);
    void readFastqReads(std::ostream &subset, std::ostream &bed, ReadTlStats &stats);

};

#endif /* INPUT_H */
