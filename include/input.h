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
    std::vector<std::string> mode;
    // uint64_t absPos = 0; // Safe for a single thread only
};


class Input {
    
    UserInputTeloscope userInput;
    std::shared_ptr<std::istream> stream;

    std::map<unsigned int, uint64_t> pathAbsPos; // Maps path ID to its absPos
    std::mutex pathAbsPosMutex; // Mutex for thread-safe access to pathAbsPos
    
public:

    std::vector<Log> logs;

    void load(UserInputTeloscope userInput);
    void read(InSequences &inSequence);
    bool walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);

};

#endif /* INPUT_H */