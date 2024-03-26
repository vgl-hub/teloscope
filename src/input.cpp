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
        // threadPool.queueJob([&inPath, this, inSegments, inGaps]() {
        threadPool.queueJob([&inPath, this, inSegments, inGaps, &teloscope]() {
            return teloscope.walkPath(&inPath, *inSegments, *inGaps); // add counter to track job order submission to be collected next, create shared object e.g. std:vector accessible to all threads. How? 
        }); 

    std::cout << "Waiting for jobs to complete" << std::endl;
    jobWait(threadPool); // Wait for all jobs to complete
    std::cout << "All jobs completed" << std::endl;
    
    // teloscope.sortWindowsBySeqPos();
    teloscope.printAllWindows();

    // Giulio: create an instance of the object to store the output. At the end of the loop, add the object to the shared vector.
    // Jack: We need a container for storing and a method for sorting the objects.
    // protect the push_back?      std::unique_lock<std::mutex> lck (mtx); std::lock_guard<std::mutex> guard(pathAbsPosMutex);
    // here put the push_back to the vector, or the shared object, to keep track of the order of jobs submitted.
}


bool Teloscope::walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps) {

    // unsigned int cUId = 0, gapLen = 0, seqPos = path->getSeqPos();
    unsigned int cUId = 0, gapLen = 0;
    std::vector<PathComponent> pathComponents = path->getComponents();
    uint64_t absPos = 0;
    std::vector<WindowData> pathWindows;

    for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
        cUId = component->id;
    
        if (component->componentType == SEGMENT) {
            
            auto inSegment = find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
            std::string sequence = (*inSegment)->getInSequence(component->start, component->end);
            unmaskSequence(sequence);
            
            if (component->orientation == '+') {
                // Jack: header = cleanString(path->getHeader()); 
                // std::vector<WindowData> segmentWindows = teloscope.analyzeSegment(path->getHeader(), sequence, userInput, absPos, pathId);
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
    insertWindowData(pathWindows);
    // insertWindowData(seqPos, pathWindows);

    return true;
}