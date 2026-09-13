#ifndef INPUT_H
#define INPUT_H

#include <stdint.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>

#include "log.h"
#include "struct.h"

class InSequences;

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
    uint8_t editDistance = 1;
    uint8_t kmerLen = 21;

    uint32_t maxMatchDist = 50;
    uint32_t minBlockLen = 300;
    bool minBlockLenSet = false;
    uint32_t maxBlockDist = 500;
    uint32_t minBlockCounts = 2;
    uint32_t terminalTolerance = 3000;
    uint32_t linkDistance = 1000;
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
    bool fastqSubset = false;
    bool bamSubset = false;

    double maxMem = 0;
    std::string prefix = ".", outFile = "";
};

bool isGfaAssemblyPath(const std::string &path);

class Input {

    UserInputTeloscope userInput;
    std::shared_ptr<std::istream> stream;
    
public:

    std::vector<Log> logs;

    void load(UserInputTeloscope userInput);
    void read(InSequences &inSequence);
    void readFastqSubset(std::ostream &out);

};

#endif /* INPUT_H */
