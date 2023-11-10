#ifndef INPUT_H
#define INPUT_H

// Add poly+inh of UserInputTeloscope

struct UserInputTeloscope : UserInput {

    std::vector<std::string> patterns = {"TTAGGG", "CCCCT"};
    int windowSize = 500;  
    uint8_t kmerLen = 21;
    int step = 1;    
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