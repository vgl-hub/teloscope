// giulio: error cout string, directly
// L40:42 std::vector<> for auto, * symbol not necessary

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

    for (InPath& inPath : inPaths) 
        walkPath(&inPath, *inSegments, *inGaps);

}

void Input::walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps) {
    
    unsigned int cUId = 0, gapLen = 0;
    
    std::vector<PathComponent> pathComponents = path->getComponents();

    uint64_t absPos = 0; // jack: put in line 52
    
    for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
        cUId = component->id;
    
        if (component->type == SEGMENT) {
            
            auto inSegment = find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
                
            std::string sequence = (*inSegment)->getInSequence(component->start, component->end);

            unmaskSequence(sequence);

            absPos += sequence.size();
                
            if (component->orientation == '+') {

                findTelomeres(path->getHeader(), sequence, userInput);

            }else{
                
                
                    

                
            }
            
        }else if (component->type == GAP){
            
            auto inGap = find_if(inGaps.begin(), inGaps.end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
            
            gapLen += inGap->getDist(component->start - component->end);

            absPos += gapLen;
            
        }else{} // need to handle edges, cigars etc
        
    }
    
}


// #include "threadpool.h" // Include the thread pool header

// void Input::walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps) {
    
//     unsigned int cUId = 0, gapLen = 0;
//     std::vector<PathComponent> pathComponents = path->getComponents();
//     uint64_t absPos = 0;
    
//     ThreadPool pool(4); // Example: Create a thread pool with 4 threads

//     for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
        
//         pool.enqueue([&, component]() { // Capture by reference & and the current component by value
//             cUId = component->id;
        
//             if (component->type == SEGMENT) {
                
//                 auto inSegment = find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;});
                    
//                 std::string sequence = (*inSegment)->getInSequence(component->start, component->end);

//                 unmaskSequence(sequence); // Apply the new function

//                 absPos += sequence.size();
                    
//                 if (component->orientation == '+') {
//                     findTelomeres(path->getHeader(), sequence, userInput);
//                 } else {
//                     // Process for '-' orientation
//                 }
                
//             } else if (component->type == GAP) {
                
//                 auto inGap = find_if(inGaps.begin(), inGaps.end(), [cUId](InGap& obj) {return obj.getuId() == cUId;});
                
//                 gapLen += inGap->getDist(component->start - component->end);

//                 absPos += gapLen;
                
//             } else {
//                 // Handle edges, cigars etc
//             }
//         });
//     }

//     pool.jobWait(); // Wait for all jobs in the thread pool to complete
// }
