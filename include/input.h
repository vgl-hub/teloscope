#ifndef INPUT_H
#define INPUT_H

#include "main.h" // not in Mac's code
#include <stdint.h>

struct UserInputTeloscope : UserInput {
    
    std::string outRoute;
    std::unordered_map<std::string, uint8_t> hammingDistances;

    std::string canonicalFwd = "CCCTAA";
    std::string canonicalRev = "TTAGGG";
    unsigned short int canonicalSize = 6;
    std::vector<std::string> patterns = {"TTAGGG", "CCCTAA"};

    uint32_t windowSize = 1000;
    uint8_t kmerLen = 21;
    uint32_t step = 500;
    unsigned short int minBlockLen = 500; // Only for all blocks
    unsigned short int maxBlockDist = 200;
    unsigned short int minBlockCounts = 2;
    unsigned short int maxMatchDist = 50;
    uint32_t terminalLimit = 50000;

    bool outFasta = false;
    bool outWinRepeats = false;
    bool outGC = false;
    bool outEntropy = false;
    bool outMatches = false;
    bool outITS = false;
    bool ultraFastMode = false;
    
    double maxMem = 0;
    std::string prefix = ".", outFile = ""; // JACK: CHECK
};

class Input {

    UserInputTeloscope userInput;
    std::shared_ptr<std::istream> stream;
    
public:

    std::vector<Log> logs;

    void load(UserInputTeloscope userInput);
    void read(InSequences &inSequence);

};

#endif /* INPUT_H */