#include <iostream>
#include <stdint.h>
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
    
    teloscope.sortBySeqPos();
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
    uint64_t pathSize = path->getLen();

    threadLog.add("\n\tWalking path:\t" + path->getHeader());
    std::string header = path->getHeader();
    eraseChar(header, '\r');
    uint32_t numSegments = (pathComponents.size() + 1) / 2;

    // Initialize PathData for this path
    PathData pathData;
    pathData.seqPos = seqPos;
    pathData.header = header;
    pathData.windows.reserve((pathSize - (userInput.windowSize + 2 * userInput.step) * (numSegments)) / userInput.step);

    for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
        
        cUId = component->id;
    
        if (component->componentType == SEGMENT) {
            
            auto inSegment = find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
            std::string sequence = (*inSegment)->getInSequence(component->start, component->end);
            unmaskSequence(sequence);
            
            if (component->orientation == '+') {
                SegmentData segmentData = analyzeSegment(sequence, userInput, absPos);

                // Collect window data
                pathData.windows.insert(
                    pathData.windows.end(),
                    std::make_move_iterator(segmentData.windows.begin()),
                    std::make_move_iterator(segmentData.windows.end())
                );

                // Collect blocks
                for (auto& [groupName, blocks] : segmentData.mergedBlocks) {
                    pathData.mergedBlocks[groupName].insert(
                        pathData.mergedBlocks[groupName].end(),
                        std::make_move_iterator(blocks.begin()),
                        std::make_move_iterator(blocks.end())
                    );
                }

                // Collect matches
                pathData.canonicalMatches.insert(
                    pathData.canonicalMatches.end(),
                    std::make_move_iterator(segmentData.canonicalMatches.begin()),
                    std::make_move_iterator(segmentData.canonicalMatches.end())
                );

                pathData.nonCanonicalMatches.insert(
                    pathData.nonCanonicalMatches.end(),
                    std::make_move_iterator(segmentData.nonCanonicalMatches.begin()),
                    std::make_move_iterator(segmentData.nonCanonicalMatches.end())
                );

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
    allPathData.push_back(std::move(pathData));

    threadLog.add("\tCompleted walking path:\t" + path->getHeader());
    logs.push_back(threadLog);

    return true;
}