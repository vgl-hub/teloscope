#ifndef INPUT_H
#define INPUT_H

#include "main.h" //check

struct UserInputTeloscope : UserInput {

    // std::vector<std::string> patterns = {"TTAGGG", "CCCTAA"};
    std::vector<std::string> patterns;
    uint32_t windowSize = 1000;
    uint8_t kmerLen = 21;
    uint32_t step = 500;
    double maxMem = 0;
    std::string prefix = ".", outFile = ""; // Jack: this should replace outRoute
    std::vector<std::string> mode;
};


class Input {
    
    UserInputTeloscope userInput;
    std::shared_ptr<std::istream> stream;
    // giulio: add vector to keep track of the order of jobs submitted. Avoid concurrency. 
    // Make sure the vector updates could be atomic. The operation  is atomic, e.g. push_back. 
    // std::vector<ThreadOutput> outputs;
    // std::mutex outputsMutex;
    
public: // Jack: indentation? 

    std::vector<Log> logs;
    
    void load(UserInputTeloscope userInput);
    
    void read(InSequences &inSequence);

    bool walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);

    // void sortOutput(); // Jack: gfalibs sort? 
    
};

#endif /* INPUT_H */