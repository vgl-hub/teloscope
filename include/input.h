#ifndef INPUT_H
#define INPUT_H

#include "main.h" // not in Mac's code

struct UserInputTeloscope : UserInput {

    std::vector<std::string> patterns;
    uint32_t windowSize = 1000;
    uint8_t kmerLen = 21;
    uint32_t step = 500;
    double maxMem = 0;
    std::string prefix = ".", outFile = "";
    bool keepWindowData = false; // Memory intensive
    bool modeMatch = true, modeEntropy = true, modeGC = true; // Change to: de novo, user-defined
    float densityThreshold = 0.5f; // Threshold for telomere detection
    uint32_t maxGaps = 2;   // Allowed gaps in telomere blocks
    uint8_t mergeDistance = 6;

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