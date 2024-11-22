#ifndef INPUT_H
#define INPUT_H

#include "main.h" // not in Mac's code

struct UserInputTeloscope : UserInput {
    
    std::string outRoute;
    std::pair<std::string, std::string> canonicalPatterns;
    std::vector<std::string> patterns;
    std::unordered_map<std::string, uint8_t> hammingDistances;

    uint32_t windowSize = 1000;
    uint8_t kmerLen = 21;
    uint32_t step = 500;
    unsigned short int minBlockLen = 100; // Not used anymore
    unsigned short int maxBlockDist = 200;
    unsigned short int minBlockCounts = 2;

    bool keepWindowData = false; // Memory intensive
    bool modeMatch = true, modeEntropy = true, modeGC = true; // Change to: de novo, user-defined

    bool outGC = true;
    bool outEntropy = true;
    bool outDistance = true;
    bool outFasta = true;
    bool outMatches = true;
    bool outITS = true;

    
    double maxMem = 0;
    std::string prefix = ".", outFile = "";
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