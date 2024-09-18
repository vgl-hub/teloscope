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
    lg.verbose("Finished loading genome assembly");

    std::vector<InPath> inPaths = inSequences.getInPaths(); 
    std::vector<InSegment*> *inSegments = inSequences.getInSegments(); 
    std::vector<InGap> *inGaps = inSequences.getInGaps();

    Teloscope teloscope(userInput);

    for (InPath& inPath : inPaths)
        threadPool.queueJob([&inPath, this, inSegments, inGaps, &teloscope]() {
            return teloscope.walkPath(&inPath, *inSegments, *inGaps);
        }); 

    lg.verbose("Waiting for jobs to complete");
    std::cout << "Waiting for jobs to complete" << std::endl;

    jobWait(threadPool); // Wait for all jobs to complete
    lg.verbose("All jobs completed");
    std::cout << "All jobs completed" << std::endl;
    
    teloscope.sortWindowsBySeqPos();
    lg.verbose("\nPaths sorted by original position");  

    teloscope.handleBEDFile();
    lg.verbose("\nBED/BEDgraph files generated");

    teloscope.printSummary();
    lg.verbose("\nSummary printed");
}


bool Teloscope::walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps) {

    Log threadLog;
    uint64_t absPos = 0;
    unsigned int cUId = 0, gapLen = 0, seqPos = path->getSeqPos();
    std::vector<PathComponent> pathComponents = path->getComponents();
    
    // std::string header = removeCarriageReturns(path->getHeader());
    std::string header = path->getHeader();
    std::cout << "header_before: " << header << std::endl;
    eraseChar(header, '\r');
    std::cout << "header_after: " << header << std::endl;
    threadLog.add("\n\tWalking path:\t" + path->getHeader());

    std::vector<WindowData> pathWindows;
    std::vector<TelomereBlock> pathTelomereBlocks;

    for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
        
        cUId = component->id;
    
        if (component->componentType == SEGMENT) {
            
            auto inSegment = find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
            std::string sequence = (*inSegment)->getInSequence(component->start, component->end);
            unmaskSequence(sequence);
            
            if (component->orientation == '+') {

                std::vector<WindowData> segmentWindows = analyzeSegment(sequence, userInput, absPos);
                pathWindows.insert(pathWindows.end(), segmentWindows.begin(), segmentWindows.end());

                // std::vector<TelomereBlock> segmentTelomereBlocks = analyzeSegment(sequence, userInput, absPos);
                // pathTelomereBlocks.insert(pathTelomereBlocks.end(), segmentTelomereBlocks.begin(), segmentTelomereBlocks.end());


            } else {
            }
            
            absPos += sequence.size();
            
        }else if (component->componentType == GAP){
            
            auto inGap = find_if(inGaps.begin(), inGaps.end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
            gapLen = inGap->getDist(component->start - component->end);

            absPos += gapLen;
            
        } else {
        } // need to handle edges, cigars etc
        
    }

    std::lock_guard<std::mutex> lck(mtx);
    insertWindowData(seqPos, header, pathWindows);

    allTelomereBlocks.push_back({seqPos, header, pathTelomereBlocks});

    // uint8_t mergeDistance = userInput.mergeDistance;
    // mergeTelomereBlocks(pathTelomereBlocks, mergeDistance);


    threadLog.add("\tCompleted walking path:\t" + path->getHeader());
    logs.push_back(threadLog);

    return true;
}