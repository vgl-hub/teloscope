#ifndef INPUT_H
#define INPUT_H

#include "main.h"

// Add poly+inh of UserInputTeloscope

struct UserInputTeloscope : UserInput { // jack: why do we need inheritance?

    std::vector<std::string> patterns = {"TTAGGG", "CCCTAA"};
    uint32_t windowSize = 500;  
    uint8_t kmerLen = 21;
    uint32_t step = 200;    
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

    void walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);
    
};

#endif /* INPUT_H */