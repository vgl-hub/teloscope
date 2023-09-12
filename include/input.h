#ifndef INPUT_H
#define INPUT_H

class Input {
    
    UserInput userInput;
    
    std::shared_ptr<std::istream> stream;
    
public:

    std::vector<Log> logs;
    
    void load(UserInput userInput);
    
    void read(InSequences &inSequence);

    void walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);
    
};

#endif /* INPUT_H */