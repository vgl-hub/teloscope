#include <iostream>
#include <algorithm>
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

namespace {

struct PendingTelomereAnnotation {
    std::string dedupeKey;
    std::string header;
    bool atStart = false;
    uint32_t blockLen = 0;
    char segOrient = '+';
};

std::vector<Tag> makeSyntheticTelomereTags(uint32_t blockLen) {
    return {
        Tag{'i', "LN", "6"},
        Tag{'i', "RC", "6000"},
        Tag{'i', "TL", std::to_string(blockLen)}
    };
}

void queueAnnotation(std::vector<PendingTelomereAnnotation> &annotations,
                     const std::string &dedupeKey,
                     const std::string &header,
                     bool atStart,
                     uint32_t blockLen,
                     char segOrient) {
    for (auto &ann : annotations) {
        if (ann.dedupeKey == dedupeKey) {
            ann.blockLen = std::max(ann.blockLen, blockLen);
            return;
        }
    }

    annotations.push_back({dedupeKey, header, atStart, blockLen, segOrient});
}

void appendTelomereConnection(InSequences &inSequences,
                              unsigned int telomereUid,
                              unsigned int segmentUid,
                              char segmentOrient) {
    InGap jump;
    jump.newGap(0, telomereUid, segmentUid, '+', segmentOrient, 0, "",
                std::vector<Tag>{Tag{'i', "RC", "0"}});
    std::lock_guard<std::mutex> lck(mtx);
    inSequences.addGap(jump);
}

} // namespace


void Input::load(UserInputTeloscope userInput) {
    
    this->userInput = userInput;
    
}


void Input::read(InSequences &inSequences) {
    loadGenome(userInput, inSequences);
    lg.verbose("Finished loading genome assembly");

    std::vector<InPath> inPaths = inSequences.getInPaths();
    std::vector<InSegment*> *inSegments = inSequences.getInSegments();
    std::vector<InGap> *inGaps = inSequences.getInGaps();
    Teloscope teloscope(userInput);

    // GFA-based annotation
    if (userInput.inSequence.find(".gfa") != std::string::npos) {

        if (!inPaths.empty()) {
            // Path-aware mode: only annotate path-terminal segments
            for (InPath& path : inPaths) {
                std::vector<PathComponent> components = path.getComponents();

                // Find first SEGMENT component
                for (const auto& comp : components) {
                    if (comp.componentType == SEGMENT && comp.orientation != '0') {
                        InSegment* seg = inSequences.getInSegment(comp.id);
                        char orient = comp.orientation;
                        threadPool.queueJob([seg, &inSequences, &teloscope, orient]() {
                            return teloscope.walkSegmentForPath(seg, inSequences, orient, true);
                        });
                        break;
                    }
                }

                // Find last SEGMENT component (reverse iteration)
                for (auto it = components.rbegin(); it != components.rend(); ++it) {
                    if (it->componentType == SEGMENT && it->orientation != '0') {
                        InSegment* seg = inSequences.getInSegment(it->id);
                        char orient = it->orientation;
                        threadPool.queueJob([seg, &inSequences, &teloscope, orient]() {
                            return teloscope.walkSegmentForPath(seg, inSequences, orient, false);
                        });
                        break;
                    }
                }
            }
        } else {
            // Fallback: path-less GFA, scan all segments (implicit + orientation)
            for (InSegment* inSegment : *inSegments) {
                threadPool.queueJob([inSegment, &inSequences, &teloscope]() {
                    return teloscope.walkSegment(inSegment, inSequences);
                });
            }
        }

        lg.verbose("Waiting for telomere annotation jobs to complete");
        jobWait(threadPool);
        lg.verbose("\nAll telomere annotation jobs completed.");

        // Create L lines from path adjacency for BandageNG graph layout
        for (InPath& path : inSequences.getInPaths()) {
            std::vector<PathComponent> components = path.getComponents();
            unsigned int prevSegId = 0;
            char prevOrient = '+';
            bool hasPrev = false;

            for (const auto& comp : components) {
                if (comp.componentType == SEGMENT && comp.orientation != '0') {
                    if (hasPrev) {
                        InEdge edge;
                        edge.newEdge(0, prevSegId, comp.id,
                                     prevOrient, comp.orientation,
                                     "0M", "", {});
                        inSequences.appendEdge(edge);
                    }
                    prevSegId = comp.id;
                    prevOrient = comp.orientation;
                    hasPrev = true;
                }
            }
        }

        // Write annotated GFA
        Report report;
        std::string outGfa = userInput.outRoute + "/" + userInput.inSequenceName + ".telo.annotated.gfa";
        report.writeToStream(inSequences, outGfa, userInput);
        lg.verbose("Annotated GFA written to " + outGfa);
        return;
    }

    // path-based annotation
    for (InPath& inPath : inPaths) {
        InPath* pathPtr = &inPath;
        threadPool.queueJob([pathPtr, inSegments, inGaps, &teloscope]() {
            return teloscope.walkPath(pathPtr, *inSegments, *inGaps);
        });
    }
    lg.verbose("Waiting for jobs to complete");
    jobWait(threadPool);
    lg.verbose("\nAll jobs completed.");

    teloscope.sortBySeqPos();
    lg.verbose("\nPaths sorted by original position.");

    teloscope.handleBEDFile();
    lg.verbose("\nReport and BED/BEDgraph files generated.");
}


bool Teloscope::walkSegment(InSegment* segment, InSequences& inSequences) {
    Log threadLog;
    threadLog.add("\n\tWalking segment:\t" + segment->getSeqHeader());

    std::string sequence = segment->getInSequence(0, 0);
    unmaskSequence(sequence);

    SegmentData segmentData = scanSegment(sequence, 0, true); // tipsOnly = true for GFA segments

    std::vector<PendingTelomereAnnotation> annotations;

    for (const TelomereBlock& block : segmentData.terminalBlocks) {
        uint64_t distToStart = block.start;
        uint64_t distToEnd   = sequence.size() - (block.start + block.blockLen);
        bool atStart = distToStart <= distToEnd;

        // edge orientation: + = start, - = end
        char segOrient = atStart ? '+' : '-';
        const std::string header = "telomere_"
                        + segment->getSeqHeader()
                        + (atStart ? "_start" : "_end");
        queueAnnotation(annotations, header, header, atStart, block.blockLen, segOrient);
    }

    for (auto& ann : annotations) {
        Sequence teloSeq{ann.header, "", new std::string("*")};
        inSequences.traverseInSegmentWrapper(&teloSeq, makeSyntheticTelomereTags(ann.blockLen));

        unsigned int teloUid;
        {
            std::lock_guard<std::mutex> lck(mtx);
            teloUid = inSequences.getHash1()->at(ann.header);
        }

        appendTelomereConnection(inSequences, teloUid, segment->getuId(), ann.segOrient);
    }

    threadLog.add("\tCompleted walking segment:\t" + segment->getSeqHeader());
    {
        std::lock_guard<std::mutex> lck(mtx);
        logs.push_back(threadLog);
    }

    return true;
}


bool Teloscope::walkSegmentForPath(InSegment* segment, InSequences& inSequences,
                                   char pathOrient, bool isFirst) {
    Log threadLog;
    threadLog.add("\n\tWalking segment (path-aware):\t" + segment->getSeqHeader());

    std::string sequence = segment->getInSequence(0, 0);
    unmaskSequence(sequence);

    SegmentData segmentData = scanSegment(sequence, 0, true);

    // Which physical end of the segment has the path-terminal tip?
    // First+'+' and Last+'-' => scan START; First+'-' and Last+'+' => scan END
    bool scanStart = (isFirst == (pathOrient == '+'));

    std::vector<PendingTelomereAnnotation> annotations;

    for (const TelomereBlock& block : segmentData.terminalBlocks) {
        uint64_t distToStart = block.start;
        uint64_t distToEnd   = sequence.size() - (block.start + block.blockLen);
        bool blockAtStart = distToStart <= distToEnd;

        // Only keep the block at the path-terminal end
        if (blockAtStart != scanStart) continue;

        std::string posLabel = isFirst ? "start" : "end";
        std::string header = "telomere_"
                        + segment->getSeqHeader()
                        + std::string(1, pathOrient)
                        + "_" + posLabel;

        // Edge orientation rule: START = pathOrient, END = flip(pathOrient)
        char edgeOrient = isFirst ? pathOrient : (pathOrient == '+' ? '-' : '+');
        queueAnnotation(annotations, header, header, blockAtStart, block.blockLen, edgeOrient);
    }

    for (auto& ann : annotations) {
        // Deduplicate across paths: check if this telomere node already exists
        {
            std::lock_guard<std::mutex> lck(mtx);
            auto* hash = inSequences.getHash1();
            if (hash->find(ann.header) != hash->end()) continue;
        }

        Sequence teloSeq{ann.header, "", new std::string("*")};
        inSequences.traverseInSegmentWrapper(&teloSeq, makeSyntheticTelomereTags(ann.blockLen));

        unsigned int teloUid;
        {
            std::lock_guard<std::mutex> lck(mtx);
            teloUid = inSequences.getHash1()->at(ann.header);
        }

        appendTelomereConnection(inSequences, teloUid, segment->getuId(), ann.segOrient);
    }

    threadLog.add("\tCompleted walking segment (path-aware):\t" + segment->getSeqHeader());
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

    // index for O(1) lookups
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

            pathData.gapInfos.push_back({absPos, static_cast<uint32_t>(gapLen)});
            absPos += gapLen;
            
        } else {
        } // need to handle edges, cigars etc
        
    }

    // Filter blocks
    labelTerminalBlocks(pathData.terminalBlocks, static_cast<uint16_t>(pathData.gapInfos.size()),
                        pathData.terminalLabel, pathData.scaffoldType,
                        pathData.pathSize, userInput.terminalLimit);
    threadLog.add("\tCompleted walking path:\t" + path->getHeader());

    std::lock_guard<std::mutex> lck(mtx);
    allPathData.push_back(std::move(pathData));
    logs.push_back(threadLog);

    return true;
}
