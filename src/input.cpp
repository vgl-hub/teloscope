#include <iostream>
#include <stdlib.h>
#include <string>
#include <stdexcept> // jack: std::runtime_error

#include "log.h"
#include "global.h"
#include "uid-generator.h"

#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "stream-obj.h"
#include "fastx.h"

#include "gfa-lines.h"
#include "gfa.h"
#include "input-gfa.h"

#include "threadpool.h"
#include "teloscope.h"
#include "input.h"

void Input::load(UserInputTeloscope userInput) {
    
    this->userInput = userInput;
    
}


void Input::read(InSequences &inSequences) {

    loadGenome(userInput, inSequences); // load from FA/FQ/GFA to templated object

    std::vector<InPath> inPaths = inSequences.getInPaths(); 
    std::vector<InSegment*> *inSegments = inSequences.getInSegments(); 
    std::vector<InGap> *inGaps = inSequences.getInGaps();

    Teloscope teloscope(userInput);

    for (InPath& inPath : inPaths)
    
        threadPool.queueJob([&inPath, this, inSegments, inGaps, &teloscope]() {
            return teloscope.walkPath(&inPath, *inSegments, *inGaps);
        }); 

    std::cout << "Waiting for jobs to complete" << "\n";
    jobWait(threadPool); // Wait for all jobs to complete
    std::cout << "All jobs completed" << "\n";
    
    teloscope.sortWindowsBySeqPos();
    teloscope.printAllWindows();
    teloscope.generateBEDFile();
}


bool Teloscope::walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps) {

    unsigned int cUId = 0, gapLen = 0, seqPos = path->getSeqPos();
    std::vector<PathComponent> pathComponents = path->getComponents();
    uint64_t absPos = 0;
    std::vector<WindowData> pathWindows;
    std::string header = cleanString(path->getHeader());

    for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
        cUId = component->id;
    
        if (component->componentType == SEGMENT) {
            
            auto inSegment = find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
            std::string sequence = (*inSegment)->getInSequence(component->start, component->end);
            unmaskSequence(sequence);
            
            if (component->orientation == '+') {

                std::vector<WindowData> segmentWindows = analyzeSegment(sequence, userInput, absPos);
                pathWindows.insert(pathWindows.end(), segmentWindows.begin(), segmentWindows.end());

            } else {
            }
            
            absPos += sequence.size();
            
        }else if (component->componentType == GAP){
            
            auto inGap = find_if(inGaps.begin(), inGaps.end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
            gapLen += inGap->getDist(component->start - component->end);

            absPos += gapLen;
            
        } else {
        } // need to handle edges, cigars etc
        
    }

    std::unique_lock<std::mutex> lck (mtx);
    insertWindowData(seqPos, header, pathWindows); 

    return true;
}