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

    pathId = 0;

    // for (InPath& inPath : inPaths)
    //     threadPool.queueJob([&inPath, this, inSegments, inGaps]() { 
    //         return walkPath(&inPath, *inSegments, *inGaps); // add counter to track job order submission to be collected next, create shared object e.g. std:vector accessible to all threads. How? 
    //     });

    for (InPath& inPath : inPaths) {
        pathAbsPos[pathId] = 0; // Initialize absPos for this path
        pathOrder.push_back(pathId); // Track the order of this path

        threadPool.queueJob([&inPath, this, inSegments, inGaps]() { 
            return walkPath(&inPath, *inSegments, *inGaps);
        });

        pathId++; // Increment for the next path
    }

    std::cout << "Waiting for jobs to complete" << std::endl;
    jobWait(threadPool); // Wait for all jobs to complete
} 

bool Input::walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps) {
    // create an instance of the object to store the output. At the end of the loop, add the object to the shared vector. 
    unsigned int cUId = 0, gapLen = 0;
    std::vector<PathComponent> pathComponents = path->getComponents();
        
    for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
        cUId = component->id;
    
        if (component->componentType == SEGMENT) {
            
            auto inSegment = find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
            std::string sequence = (*inSegment)->getInSequence(component->start, component->end);
            unmaskSequence(sequence);
            
            if (component->orientation == '+') {

                findTelomeres(path->getHeader(), sequence, userInput);

            } else {
            }
            
            pathAbsPos[pathId] += sequence.size();
            
        }else if (component->componentType == GAP){
            
            auto inGap = find_if(inGaps.begin(), inGaps.end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
            gapLen += inGap->getDist(component->start - component->end);

            pathAbsPos[pathId] += gapLen;
            
        } else {
        } // need to handle edges, cigars etc
        
    }

    // userInput.absPos = 0; // reset the absolute position for the next segment

    // protect the push_back?      std::unique_lock<std::mutex> lck (mtx);
    // here put the push_back to the vector, or the shared object, to keep track of the order of jobs submitted.
    return true;
}