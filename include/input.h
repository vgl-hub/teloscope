#ifndef INPUT_H
#define INPUT_H

#include "main.h"

// Add poly+inh of UserInputTeloscope

struct UserInputTeloscope : UserInput { // Inheritance to reuse UserInput properties

    std::vector<std::string> patterns = {"TTAGGG", "CCCTAA"};
    uint32_t windowSize = 1000;
    uint8_t kmerLen = 21;
    uint32_t step = 500;
    double maxMem = 0;
    std::string prefix = ".", outFile = "";
    std::vector<std::string> mode;
};


class Input {
    
    UserInputTeloscope userInput;
    
    std::shared_ptr<std::istream> stream;
    
public:

    std::vector<Log> logs;
    
    void load(UserInputTeloscope userInput);
    
    void read(InSequences &inSequence);

    void walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);
    
};

#endif /* INPUT_H */