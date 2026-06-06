#include <iostream>
#include <algorithm>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <set>
#include <tuple>
#include <fstream>
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
#include "read-filter.h"

namespace {

struct PendingTelomereAnnotation {
    std::string dedupeKey;
    std::string header;
    bool atStart = false;
    uint32_t blockLen = 0;
    char segOrient = '+';
};

struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string plus;
    std::string quality;
};

struct FastqChunkResult {
    std::string output;
    uint64_t scanned = 0;
    uint64_t passed = 0;
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
    // Telomere cap is a direct adjacency: an L link with 0M overlap, not a J gap.
    InEdge edge;
    edge.newEdge(0, telomereUid, segmentUid, '+', segmentOrient, "0M", "",
                 std::vector<Tag>{Tag{'i', "RC", "0"}});
    inSequences.appendEdge(edge);
}

[[noreturn]] void fastqExitFailure() {
    threadPool.join();
    exit(EXIT_FAILURE);
}

size_t logicalLineLength(const std::string &line) {
    if (!line.empty() && line.back() == '\r') {
        return line.size() - 1;
    }
    return line.size();
}

[[noreturn]] void fastqInputError(uint64_t recordNumber, const std::string &message) {
    if (recordNumber > 0) {
        fprintf(stderr, "Error: FASTQ record %" PRIu64 ": %s.\n",
                recordNumber, message.c_str());
    } else {
        fprintf(stderr, "Error: %s.\n", message.c_str());
    }
    fastqExitFailure();
}

bool readFastqRecord(std::istream &stream, FastqRecord &record, uint64_t recordNumber) {
    // Skip blank lines (including a lone '\r') before a header so a trailing newline
    // or a stray blank line does not abort an otherwise valid stream.
    do {
        if (!std::getline(stream, record.header)) {
            return false;
        }
    } while (logicalLineLength(record.header) == 0);
    if (!std::getline(stream, record.sequence) ||
        !std::getline(stream, record.plus) ||
        !std::getline(stream, record.quality)) {
        fastqInputError(recordNumber, "truncated FASTQ record");
    }

    if (record.header.empty() || record.header.front() != '@') {
        fastqInputError(recordNumber, "expected header line starting with '@'");
    }
    if (record.plus.empty() || record.plus.front() != '+') {
        fastqInputError(recordNumber, "expected separator line starting with '+'");
    }
    if (logicalLineLength(record.sequence) != logicalLineLength(record.quality)) {
        fastqInputError(recordNumber, "sequence and quality length differ");
    }

    return true;
}

void appendFastqRecord(std::string &out, const FastqRecord &record) {
    out += record.header;
    out += '\n';
    out += record.sequence;
    out += '\n';
    out += record.plus;
    out += '\n';
    out += record.quality;
    out += '\n';
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

        // One scan job per terminal end. Resolve every job before queuing any:
        // worker threads append telomere nodes to the live segment vector, so
        // reading or iterating it while jobs run would race that growth.
        struct TeloJob { InSegment* seg; char orient; bool isFirst; bool pathAware; };
        std::vector<TeloJob> jobs;

        if (!inPaths.empty()) {
            // Unique (segment, orientation, isFirst) terminal ends: one job per node.
            std::set<std::tuple<unsigned int, char, bool>> terminalEnds;
            for (InPath& path : inPaths) {
                std::vector<PathComponent> components = path.getComponents();

                // First SEGMENT component => contig start (isFirst = true)
                for (const auto& comp : components) {
                    if (comp.componentType == SEGMENT && comp.orientation != '0') {
                        terminalEnds.emplace(comp.id, comp.orientation, true);
                        break;
                    }
                }

                // Last SEGMENT component => contig end (isFirst = false)
                for (auto it = components.rbegin(); it != components.rend(); ++it) {
                    if (it->componentType == SEGMENT && it->orientation != '0') {
                        terminalEnds.emplace(it->id, it->orientation, false);
                        break;
                    }
                }
            }

            for (const auto& end : terminalEnds)
                jobs.push_back({inSequences.getInSegment(std::get<0>(end)),
                                std::get<1>(end), std::get<2>(end), true});
        } else {
            // Fallback: path-less GFA, scan all segments (implicit + orientation)
            for (InSegment* inSegment : *inSegments)
                jobs.push_back({inSegment, '+', false, false});
        }

        // Count missing sequence now, while the segment vector is still stable.
        size_t segmentsScanned = jobs.size(), segmentsNoSeq = 0;
        for (const TeloJob& job : jobs)
            if (job.seg->getInSequencePtr() == NULL) ++segmentsNoSeq;

        // No DNA means nothing to scan: say so instead of silently annotating nothing.
        if (segmentsNoSeq > 0)
            fprintf(stderr, "Warning: %zu of %zu GFA segment(s) had no sequence (*); skipped for telomere annotation.\n",
                    segmentsNoSeq, segmentsScanned);

        for (const TeloJob& job : jobs) {
            InSegment* seg = job.seg;
            char orient = job.orient;
            bool isFirst = job.isFirst;
            if (job.pathAware)
                threadPool.queueJob([seg, &inSequences, &teloscope, orient, isFirst]() {
                    return teloscope.walkSegmentForPath(seg, inSequences, orient, isFirst);
                });
            else
                threadPool.queueJob([seg, &inSequences, &teloscope]() {
                    return teloscope.walkSegment(seg, inSequences);
                });
        }

        lg.verbose("Waiting for telomere annotation jobs to complete");
        jobWait(threadPool);
        lg.verbose("\nAll telomere annotation jobs completed.");

        // Write annotated GFA
        Report report;
        std::string outGfa = userInput.outRoute + "/" + userInput.inSequenceName + ".telo.annotated.gfa";
        report.writeToStream(inSequences, outGfa, userInput);
        lg.verbose("Annotated GFA written to " + outGfa);

        // Companion colours CSV: every telomere node green (#008000), read from the graph.
        std::string outColors = userInput.outRoute + "/" + userInput.inSequenceName + ".telo.annotated.colors.csv";
        std::ofstream colorsStream(outColors);
        if (colorsStream) {
            colorsStream << "node\tcolor\n";
            for (InSegment* seg : *inSequences.getInSegments()) {
                const std::string& name = seg->getSeqHeader();
                if (name.rfind("telomere_", 0) == 0)
                    colorsStream << name << "\t#008000\n";
            }
            lg.verbose("Telomere colours CSV written to " + outColors);
        } else {
            fprintf(stderr, "Warning: could not write telomere colours CSV to %s.\n", outColors.c_str());
        }
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


void Input::readFastqSubset(std::ostream &out) {
    StreamObj streamObj;
    std::shared_ptr<std::istream> stream = streamObj.openStream(userInput, 'f');
    if (!stream) {
        fprintf(stderr, "Error: Stream not successful: %s.\n", userInput.inSequence.c_str());
        fastqExitFailure();
    }

    int first = stream->peek();
    if (first == EOF) {
        fastqInputError(0, "FASTQ input is empty");
    }
    if (first != '@') {
        fastqInputError(0, "FASTQ input must start with '@'");
    }

    const uint32_t threads = std::max<uint32_t>(1, threadPool.totalThreads());
    const size_t recordsPerBatch = std::min<size_t>(
        2048, std::max<size_t>(256, static_cast<size_t>(threads) * 32));

    std::vector<FastqRecord> batch;
    batch.reserve(recordsPerBatch);

    uint64_t recordNumber = 0;
    uint64_t totalReads = 0;
    uint64_t passedReads = 0;

    auto processBatch = [&]() {
        if (batch.empty()) {
            return;
        }

        const size_t chunkCount = std::min<size_t>(threads, batch.size());
        const size_t chunkSize = (batch.size() + chunkCount - 1) / chunkCount;
        std::vector<FastqChunkResult> results(chunkCount);

        for (size_t chunk = 0; chunk < chunkCount; ++chunk) {
            const size_t start = chunk * chunkSize;
            const size_t end = std::min(batch.size(), start + chunkSize);
            if (start >= end) {
                continue;
            }

            threadPool.queueJob([&, chunk, start, end]() {
                ReadTelomereFilter filter(userInput);
                FastqChunkResult &result = results[chunk];
                result.scanned = end - start;

                for (size_t i = start; i < end; ++i) {
                    if (filter.matches(batch[i].sequence)) {
                        appendFastqRecord(result.output, batch[i]);
                        result.passed++;
                    }
                }

                return true;
            });
        }

        jobWait(threadPool);

        for (const auto &result : results) {
            if (!result.output.empty()) {
                out.write(result.output.data(), static_cast<std::streamsize>(result.output.size()));
            }
            totalReads += result.scanned;
            passedReads += result.passed;
        }

        if (!out.good()) {
            fprintf(stderr, "Error: failed while writing FASTQ subset to stdout.\n");
            fastqExitFailure();
        }

        batch.clear();
    };

    while (true) {
        FastqRecord record;
        if (!readFastqRecord(*stream, record, recordNumber + 1)) {
            break;
        }
        recordNumber++;
        batch.push_back(std::move(record));

        if (batch.size() == recordsPerBatch) {
            processBatch();
        }
    }

    processBatch();
    out.flush();

    fprintf(stderr, "FASTQ subset: kept %" PRIu64 " of %" PRIu64 " reads.\n",
            passedReads, totalReads);
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
        // Path-less segments are scanned as '+', so encode '+' in the node name.
        const std::string header = "telomere_"
                        + segment->getSeqHeader()
                        + "+"
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

    // Physical end holding the path-terminal tip (start when isFirst == (orient=='+')).
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
        // Terminal ends were de-duplicated upstream, so each header is created once.
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
