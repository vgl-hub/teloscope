#ifndef INPUT_H
#define INPUT_H

#include "main.h" //check

struct UserInputTeloscope : UserInput {

    std::vector<std::string> patterns;
    uint32_t windowSize = 1000;
    uint8_t kmerLen = 21;
    uint32_t step = 500;
    double maxMem = 0;
    std::string prefix = ".", outFile = "";
    std::vector<std::string> mode;
    uint64_t absPos = 0;
};


class Input {
    
    UserInputTeloscope userInput;
    std::shared_ptr<std::istream> stream;
    
public:

    std::vector<Log> logs;

    void load(UserInputTeloscope userInput);
    void read(InSequences &inSequence);
    bool walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);

};

#endif /* INPUT_H */