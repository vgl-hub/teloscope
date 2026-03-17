#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <unordered_map>
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

    SegmentData segmentData = scanSegment(sequence, 0, true); // tipsOnly = true for GFA segments

    // Collect annotations before locking to minimize critical section
    struct TeloAnnotation {
        std::string header;
        std::vector<Tag> tags;
        char segOrient; // orientation of the assembly segment in the edge
    };
    std::vector<TeloAnnotation> annotations;

    for (const TelomereBlock& block : segmentData.terminalBlocks) {
        uint64_t distToStart = block.start;
        uint64_t distToEnd   = sequence.size() - (block.start + block.blockLen);
        bool atStart = distToStart <= distToEnd;

        std::string header = "telomere_"
                        + segment->getSeqHeader()
                        + (atStart ? "_start" : "_end");

        // Skip duplicate telomere at the same end (e.g. overlapping fwd+rev blocks)
        bool duplicate = false;
        for (const auto& ann : annotations) {
            if (ann.header == header) { duplicate = true; break; }
        }
        if (duplicate) continue;

        // Edge orientation for BandageNG:
        //   L telo + segment + 0M  →  telo at left (start) of segment
        //   L telo + segment - 0M  →  telo at right (end) of segment
        char segOrient = atStart ? '+' : '-';

        annotations.push_back({header, {
            Tag{'i', "LN", "6"},
            Tag{'i', "RC", "6000"},
            Tag{'i', "TL", std::to_string(block.blockLen)}
        }, segOrient});
    }

    // Note: traverseInSegmentWrapper and appendEdge lock the global mtx internally.
    // We must NOT hold mtx while calling them, or we deadlock.
    for (auto& ann : annotations) {
        // Create telomere segment (locks mtx internally)
        Sequence teloSeq{ann.header, "", new std::string("*")};
        inSequences.traverseInSegmentWrapper(&teloSeq, ann.tags);

        // Read back UID under lock (hash map is not thread-safe for concurrent R/W)
        unsigned int teloUid;
        {
            std::lock_guard<std::mutex> lck(mtx);
            teloUid = inSequences.getHash1()->at(ann.header);
        }

        // Create and append edge (locks mtx internally)
        InEdge inEdge;
        inEdge.newEdge(
          0,                      // auto-assign UID in appendEdge
          teloUid,                // source = telomere node
          segment->getuId(),      // target = assembly segment
          '+',                    // telomere always + orientation
          ann.segOrient,          // segment orient: + for start, - for end
          "0M",
          "",
          std::vector<Tag>{ Tag{'i',"RC","0"} }
        );
        inSequences.appendEdge(inEdge);
    }

    threadLog.add("\tCompleted walking segment:\t" + segment->getSeqHeader());
    {
        std::lock_guard<std::mutex> lck(mtx);
        logs.push_back(threadLog);
    }

    return true;
}


bool Teloscope::walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps) {
    Log threadLog;
    uint64_t absPos = 0;
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

    // Build index for O(1) lookups instead of O(N) find_if per component
    std::unordered_map<unsigned int, InSegment*> segmentIndex;
    segmentIndex.reserve(inSegments.size());
    for (auto* seg : inSegments) segmentIndex[seg->getuId()] = seg;

    std::unordered_map<unsigned int, InGap*> gapIndex;
    gapIndex.reserve(inGaps.size());
    for (auto& gap : inGaps) gapIndex[gap.getuId()] = &gap;

    for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
        cUId = component->id;

        if (component->componentType == SEGMENT) {
            auto inSegment = segmentIndex.find(cUId)->second;
            std::string sequence = inSegment->getInSequence(component->start, component->end);
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
            
            auto inGap = gapIndex.find(cUId)->second;
            gapLen = inGap->getDist(component->start - component->end);

            absPos += gapLen;
            pathData.gaps++;
            
        } else {
        } // need to handle edges, cigars etc
        
    }

    // Filter blocks
    labelTerminalBlocks(pathData.terminalBlocks, pathData.gaps,
                        pathData.terminalLabel, pathData.scaffoldType,
                        pathData.pathSize, userInput.terminalLimit);
    threadLog.add("\tCompleted walking path:\t" + path->getHeader());

    std::lock_guard<std::mutex> lck(mtx);
    allPathData.push_back(std::move(pathData));
    logs.push_back(threadLog);

    return true;
}