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
#include "output.h"
#include "input.h"


void Input::load(UserInputTeloscope userInput) {
    
    this->userInput = userInput;
    
}


void Input::read(InSequences &inSequences) {
    // Load sequences (FA, FQ, or GFA) into inSequences and prepare handlers
    loadGenome(userInput, inSequences);
    lg.verbose("Finished loading genome assembly");

    std::vector<InPath> inPaths = inSequences.getInPaths();
    std::vector<InSegment*> *inSegments = inSequences.getInSegments();
    std::vector<InGap> *inGaps = inSequences.getInGaps();
    Teloscope teloscope(userInput);

    // GFA-based annotation
    if (userInput.inSequence.find(".gfa") != std::string::npos) {
        for (InSegment* inSegment : *inSegments) {
            threadPool.queueJob([inSegment, &inSequences, &teloscope]() {
                return teloscope.walkSegment(inSegment, inSequences);
            });
        }
        
        lg.verbose("Waiting for telomere annotation jobs to complete");
        jobWait(threadPool);
        lg.verbose("\nAll telomere annotation jobs completed.");

        // Write annotated GFA
        Report report;
        std::string outGfa = userInput.outRoute + "/assembly.telo.annotated.gfa";
        report.writeToStream(inSequences, outGfa, userInput);
        lg.verbose("Annotated GFA written to " + outGfa);
        return;
    }

    // Path-based annotation and BED export
    for (InPath& inPath : inPaths) {
        threadPool.queueJob([&inPath, this, inSegments, inGaps, &teloscope]() {
            return teloscope.walkPath(&inPath, *inSegments, *inGaps);
        });
    }
    lg.verbose("Waiting for jobs to complete");
    jobWait(threadPool);
    lg.verbose("\nAll jobs completed.");

    teloscope.sortBySeqPos();
    lg.verbose("\nPaths sorted by original position.");

    teloscope.handleBEDFile();
    lg.verbose("\nBED/BEDgraph files generated.");

    teloscope.printSummary();
    lg.verbose("\nSummary printed.");
}


bool Teloscope::walkSegment(InSegment* segment, InSequences& inSequences) {
    Log threadLog;
    threadLog.add("\n\tWalking segment:\t" + segment->getSeqHeader());

    std::string sequence = segment->getInSequence(0, 0);
    unmaskSequence(sequence);

    // Initialize SegmentData for this segment
    SegmentData segmentData;
    segmentData = scanSegment(sequence, 0, true); // tipsOnly = true for GFA segments
    // segmentData.terminalBlocks = filterTerminalBlocks(segmentData.terminalBlocks); // TODO: extendBlocks???
    // char segOr = (block.blockLabel == 'p' ? '+' : '-'); // TODO
    // auto &adj = inSequences.getAdjEdgeList()[ segment->getuId() ];


    std::cout << "Sequence: " << sequence.substr(0,1000) << std::endl;
    std::cout << "Terminal blocks: " << segmentData.terminalBlocks.size() << std::endl;

    for (const TelomereBlock& block : segmentData.terminalBlocks) {
        uint64_t distToStart = block.start;
        uint64_t distToEnd   = sequence.size() - (block.start + block.blockLen);
        bool atStart = distToStart <= distToEnd;
        char segOr = (block.blockLabel == 'p' ? '+' : '-'); 

        // Build a unique header for new telomere node
        std::string header = "telomere_"
                        + segment->getSeqHeader()
                        + segOr
                        + (atStart ? "_start" : "_end");
        std::cout << "Header: " << header << std::endl;

        // Placeholder node with tags: LN (unit length), RC (read count), TL (telomere length)
        Sequence* teloSeq = new Sequence{ header, "", new std::string("*") }; // Empty sequence "*"
        std::vector<Tag> tags = {
            Tag{'i', "LN", "6"},
            Tag{'i', "RC", "6000"},
            Tag{'i', "TL", std::to_string(block.blockLen)}
        };
        std::cout << "Tags: " << tags[0].content << ", " << tags[1].content << ", " << tags[2].content << std::endl;

        // Create & hash the new telomere segment synchronously
        inSequences.traverseInSegmentWrapper(teloSeq, tags);
        unsigned int teloUid = inSequences.getHash1()->at(header);
        printf("Telomere node UID: %u\n", teloUid);

        // Prepare edge orientation
        char fromOrient = '+';
        char toOrient   = atStart
                            ? segOr
                            : (segOr == '+' ? '-' : '+');

        // Link telomere node back to our segment
        //    L telomere…      +   <seg>   <segOr>   0M   RC:i:0
        InEdge inEdge;
        inEdge.newEdge(
          inSequences.uId.next(),   // new edge UID
          segment->getuId(),            // source = our contig
          teloUid,                  // target = telomere node
          fromOrient,                      // FromOrient
          toOrient,                      // ToOrient
          "0M",                      // no overlap
          "",                       // no header
          std::vector<Tag>{ Tag{'i',"RC","0"} } // no read counts
        );
        inSequences.insertHash(inEdge.geteHeader(), inEdge.geteUId());
        inSequences.appendEdge(inEdge);
    }

    threadLog.add("\tCompleted walking segment:\t" + segment->getSeqHeader());
    logs.push_back(threadLog);
    return true;
}


bool Teloscope::walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps) {
    Log threadLog;
    uint32_t absPos = 0;
    unsigned int cUId = 0, gapLen = 0, seqPos = path->getSeqPos();
    std::vector<PathComponent> pathComponents = path->getComponents();

    threadLog.add("\n\tWalking path:\t" + path->getHeader());
    std::string header = path->getHeader();
    eraseChar(header, '\r');

    // Initialize PathData for this path
    PathData pathData;
    pathData.seqPos = seqPos;
    pathData.header = header;
    pathData.pathSize = path->getLen();
    // pathData.windows.reserve(inSegments.size()); NumWindows = ceil((L - W) / S) + 1

    for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
        cUId = component->id;
    
        if (component->componentType == SEGMENT) {
            auto inSegment = find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
            std::string sequence = (*inSegment)->getInSequence(component->start, component->end);
            unmaskSequence(sequence);
            
            if (component->orientation == '+') {
                SegmentData segmentData = scanSegment(sequence, absPos, userInput.ultraFastMode);

                // Collect window data
                pathData.windows.insert(
                    pathData.windows.end(),
                    std::make_move_iterator(segmentData.windows.begin()),
                    std::make_move_iterator(segmentData.windows.end())
                );

                // Collect blocks
                pathData.terminalBlocks.insert(
                    pathData.terminalBlocks.end(),
                    std::make_move_iterator(segmentData.terminalBlocks.begin()),
                    std::make_move_iterator(segmentData.terminalBlocks.end())
                );

                pathData.interstitialBlocks.insert(
                    pathData.interstitialBlocks.end(),
                    std::make_move_iterator(segmentData.interstitialBlocks.begin()),
                    std::make_move_iterator(segmentData.interstitialBlocks.end())
                );

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
            pathData.gaps++;
            
        } else {
        } // need to handle edges, cigars etc
        
    }

    // Filter blocks
    labelTerminalBlocks(pathData.terminalBlocks, pathData.gaps, 
                        pathData.terminalLabel, pathData.scaffoldType);
    pathData.interstitialBlocks = filterITSBlocks(pathData.interstitialBlocks);

    std::lock_guard<std::mutex> lck(mtx);
    allPathData.push_back(std::move(pathData));

    threadLog.add("\tCompleted walking path:\t" + path->getHeader());
    logs.push_back(threadLog);

    return true;
}