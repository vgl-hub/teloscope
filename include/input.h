#ifndef INPUT_H
#define INPUT_H

#include "main.h" // not in Mac's code

struct UserInputTeloscope : UserInput {
    
    std::string outRoute;
    std::pair<std::string, std::string> canonicalPatterns;
    std::vector<std::string> patterns;
    uint32_t windowSize = 1000;
    uint8_t kmerLen = 21;
    uint32_t step = 500;
    unsigned short int minBlockLen = 24;
    unsigned short int minBlockDist = 5;
    unsigned short int maxBlockDist = 100;

    bool keepWindowData = false; // Memory intensive
    bool modeMatch = true, modeEntropy = true, modeGC = true; // Change to: de novo, user-defined
    
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