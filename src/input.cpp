// Tasks: jack
// throw std::runtime_error("Invalid input");
// throw std::runtime_error("Stream not initialized");
// throw std::runtime_error("Null path");
// throw std::runtime_error("Segment not found");
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

#include "teloscope.h"
#include "input.h"

void Input::load(UserInput userInput) {
    
    this->userInput = userInput;
    
}

void Input::read(InSequences &inSequences) {

    loadSequences(userInput, &inSequences, 'f', 0); // load from FASTA/FASTQ to templated object

    jobWait(threadPool);

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

// jack attempt: simpler range-based for loop, ./-> object/pointer
void Input::walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps) {
    
    unsigned int cUId = 0, gapLen = 0;
    std::vector<PathComponent> pathComponents = path->getComponents(); // auto? 
    uint64_t absPos = 0;
    
    for (const auto& component : pathComponents) { // giulio: it doesnt knonw where you are in the vector
        cUId = component.id; // jack: not sure if this will work
    
        if (component.type == SEGMENT) {
            
            auto inSegment = std::find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment* obj) { // giulio: lambda: [ciud] space of the function
                return obj->getuId() == cUId; // giulio: looking at the obj in the genome
            });

            std::string sequence = (*inSegment)->getInSequence(component.start, component.end);
            absPos += sequence.size();

            if (component.orientation == '+') {
                findTelomeres(path->getHeader(), sequence, userInput); // giulio: what is your purpose?
            } else {} // giulio: gfa support (eventually!)
            
        } else if (component.type == GAP) {
            auto inGap = std::find_if(inGaps.begin(), inGaps.end(), [cUId](InGap& obj) {
                return obj.getuId() == cUId;
            });
            
            gapLen += inGap->getDist(component.start - component.end);
            absPos += gapLen;
        } else {}
            // need to handle edges, cigars etc //  giulio: SEGMENT, GAP, EDGE (not a priority right now)    
    }
}
