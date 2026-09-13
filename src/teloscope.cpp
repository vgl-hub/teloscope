#include <iostream>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>
#include <chrono>

#include "log.h"
#include "global.h"
#include "uid-generator.h"
#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "gfa-lines.h"
#include "gfa.h"
#include "sak.h"
#include "stream-obj.h"
#include "input-agp.h"
#include "input-filters.h"
#include "input-gfa.h"
#include "teloscope.h"
#include "input.h"


#ifndef TELOSCOPE_COMMIT
#define TELOSCOPE_COMMIT "unknown"
#endif

// per-base coverage runs; defined at file scope so include/teloscope.h can forward-declare it
struct CoverRun {
    uint64_t start;
    uint32_t len;
};

namespace {

constexpr std::string_view teloscopeVersion = "0.1.6";

void writeProvenanceHeader(std::ofstream& file, const UserInputTeloscope& input,
                           std::string_view columns) {
    if (!file.is_open()) return;

    std::ostringstream header;
    header << std::boolalpha
           << "#teloscope version=" << teloscopeVersion
           << " commit=" << TELOSCOPE_COMMIT << '\n'
           << "#params canonical=" << input.canonicalFwd << '/' << input.canonicalRev
           << " patterns=" << input.patterns.size()
           << " window=" << input.windowSize
           << " step=" << input.step
           << " terminal_limit=" << input.terminalLimit
           << " max_match_dist=" << input.maxMatchDist
           << " max_block_dist=" << input.maxBlockDist
           << " min_block_len=" << input.minBlockLen
           << " min_block_density=" << input.minBlockDensity
           << " min_block_counts=" << input.minBlockCounts
           << " min_canonical_count=4"
           << " terminal_tolerance=" << input.terminalTolerance
           << " link_distance=" << input.linkDistance
           << " label_threshold=" << input.labelThreshold
           << " edit_distance=" << static_cast<unsigned int>(input.editDistance)
           << " ultra_fast=" << input.ultraFastMode
           << " manual_curation=" << input.manualCuration << '\n'
           << "#columns\t" << columns << '\n';
    file << header.str();
}

const char* junctionToString(char code) {
    switch (code) {
        case 'f': return "fusion";
        case 't': return "tail_to_tail";
        case 'r': return "fragmentation";
        default:  return "single";
    }
}

char closestEnd(uint64_t start, uint32_t blockLen, uint64_t pathSize, char anchorSide) {
    if (anchorSide != '\0') return anchorSide; // terminal row: the end its telomere belongs to
    uint64_t mid2 = 2 * start + blockLen;
    return (mid2 <= pathSize) ? 'p' : 'q';
}

// R9: chr start end teloLen teloLabel closestEnd fwdCan revCan fwdNonCan revNonCan chrSize teloType
void writeBlockRow(std::ofstream& file, std::string_view pathName,
                   const TelomereBlock& block, uint64_t pathSize, std::string_view teloType) {
    file << pathName << '\t'
         << block.start << '\t'
         << (block.start + block.blockLen) << '\t'
         << block.teloLen << '\t'
         << block.strandLabel << '\t'
         << closestEnd(block.start, block.blockLen, pathSize, block.anchorSide) << '\t'
         << block.fwdCanCount << '\t'
         << block.revCanCount << '\t'
         << block.fwdNonCanCount << '\t'
         << block.revNonCanCount << '\t'
         << pathSize << '\t'
         << teloType << '\n';
}

constexpr uint32_t minCanonicalCount = 4;

enum class Orient { Fwd, Rev, All };

struct Seed {
    uint64_t start;
    uint64_t end;
    uint64_t counts;
    bool isForward;
};

// per-base coverage, overlapping matches OR-ed instead of summed
void getCoverRuns(const std::vector<MatchInfo>& matches, Orient orient, std::vector<CoverRun>& runs) {
    for (const MatchInfo& match : matches) {
        bool keep;
        switch (orient) {
            case Orient::Fwd: keep = match.isCanonical && match.isForward; break;
            case Orient::Rev: keep = match.isCanonical && !match.isForward; break;
            default:          keep = true; break;
        }
        if (!keep) continue;
        uint64_t end = match.position + match.matchSize;
        if (!runs.empty() && match.position <= runs.back().start + runs.back().len) {
            if (end > runs.back().start + runs.back().len) {
                runs.back().len = static_cast<uint32_t>(end - runs.back().start);
            }
        } else {
            runs.push_back({match.position, static_cast<uint32_t>(match.matchSize)});
        }
    }
}

// a lone opposite-strand hit is noise inside an array, a run of them is a junction
constexpr uint32_t inversionRun = 3;

void getSeeds(std::vector<MatchInfo>::const_iterator first, std::vector<MatchInfo>::const_iterator last,
              uint32_t matchDist, std::vector<Seed>& seeds) {
    bool seedFwd = false;
    uint64_t oppStart = 0, oppEnd = 0, oppCounts = 0, keptEnd = 0;
    uint32_t oppLen = 0;

    for (auto it = first; it != last; ++it) {
        const MatchInfo& match = *it;
        uint64_t end = match.position + match.matchSize;

        if (seeds.empty() || match.position > seeds.back().end + matchDist) {
            seeds.push_back({match.position, end, 1, static_cast<bool>(match.isForward)});
            seedFwd = match.isForward;
            oppLen = 0;
            continue;
        }

        if (static_cast<bool>(match.isForward) != seedFwd) {
            if (oppLen == 0) {
                oppStart = match.position;
                oppEnd = end;
                oppCounts = 0;
                keptEnd = seeds.back().end;
            }
            ++oppLen;
            ++oppCounts;
            if (end > oppEnd) oppEnd = end;
        } else {
            oppLen = 0;
        }

        if (end > seeds.back().end) seeds.back().end = end;
        seeds.back().counts++;

        if (oppLen == inversionRun) { // hand the inverted run to a seed of its own
            // the cut has to leave the seeds disjoint, or the junction clamp misses them
            seeds.back().end = std::min(keptEnd, oppStart);
            seeds.back().counts -= oppCounts;
            seeds.push_back({oppStart, oppEnd, oppCounts, static_cast<bool>(match.isForward)});
            seedFwd = match.isForward;
            oppLen = 0;
        }
    }
}
// runs are ascending, so a caller walking ascending ranges keeps one cursor
uint64_t getCoveredBases(const std::vector<CoverRun>& runs, uint64_t from, uint64_t to, size_t& cursor) {
    uint64_t covered = 0;
    while (cursor < runs.size() && runs[cursor].start + runs[cursor].len <= from) ++cursor;
    for (size_t i = cursor; i < runs.size() && runs[i].start < to; ++i) {
        uint64_t start = std::max(from, runs[i].start), end = std::min(to, runs[i].start + runs[i].len);
        if (end > start) covered += end - start;
    }
    return covered;
}

void tallyBlock(const std::vector<MatchInfo>& matches, TelomereBlock& block, size_t& cursor) {
    uint64_t end = block.start + block.blockLen;
    while (cursor < matches.size() && matches[cursor].position < block.start) ++cursor;
    for (size_t i = cursor; i < matches.size() && matches[i].position < end; ++i) {
        if (matches[i].position < block.start) continue; // cursor is shared, order is not assumed
        if (matches[i].isCanonical) {
            if (matches[i].isForward) block.fwdCanCount++; else block.revCanCount++;
        } else {
            if (matches[i].isForward) block.fwdNonCanCount++; else block.revNonCanCount++;
        }
    }
}

// one-off caller: seek with a binary search instead of a cursor started at 0
void tallyBlock(const std::vector<MatchInfo>& matches, TelomereBlock& block) {
    uint64_t end = block.start + block.blockLen;
    auto it = std::lower_bound(matches.begin(), matches.end(), block.start,
        [](const MatchInfo& m, uint64_t v) { return m.position < v; });
    for (; it != matches.end() && it->position < end; ++it) {
        if (it->isCanonical) {
            if (it->isForward) block.fwdCanCount++; else block.revCanCount++;
        } else {
            if (it->isForward) block.fwdNonCanCount++; else block.revNonCanCount++;
        }
    }
}

// walk inward from the anchor and cut where the cumulative score peaks (no gap bookkeeping: a contig has none)
uint64_t trimInward(const std::vector<CoverRun>& runs, uint64_t anchor, uint64_t limit,
                    uint32_t maxBlockDist, float weight, bool toRight, uint64_t& probe) {
    double cumulative = 0.0;
    uint64_t bestPos = anchor, prev = anchor;
    bool bridging = false; // only true once a covered run is behind us
    probe = anchor; // furthest run edge visited (R6 probe); unchanged if none

    if (toRight) {
        for (const CoverRun& run : runs) {
            if (run.start + run.len <= anchor) continue;
            probe = std::max(probe, run.start + run.len);
            if (run.start >= limit) break;
            uint64_t from = std::max(run.start, anchor);
            if (bridging && (from - prev) > maxBlockDist) break;
            uint64_t end = std::min<uint64_t>(run.start + run.len, limit);
            cumulative -= weight * static_cast<double>(from - prev);
            if (end > from) cumulative += end - from;
            if (cumulative >= 0.0) bestPos = end;
            bridging = true;
            prev = run.start + run.len;
        }
    } else {
        for (auto run = runs.rbegin(); run != runs.rend(); ++run) {
            uint64_t end = run->start + run->len;
            if (run->start >= anchor) continue;
            probe = std::min(probe, run->start);
            if (end <= limit) break;
            uint64_t to = std::min(end, anchor);
            if (bridging && (prev - to) > maxBlockDist) break;
            uint64_t start = std::max(run->start, limit);
            cumulative -= weight * static_cast<double>(prev - to);
            if (to > start) cumulative += to - start;
            if (cumulative >= 0.0) bestPos = start;
            bridging = true;
            prev = run->start;
        }
    }
    return bestPos;
}

// every locally maximal qualifying segment, so neighbouring arrays stay separate
void getMaximalSegments(const std::vector<CoverRun>& runs, uint64_t from, uint64_t to,
                        uint32_t maxBlockDist, float weight, size_t& cursor,
                        std::vector<std::pair<uint64_t, uint64_t>>& segments) {
    double cumulative = 0.0;
    uint64_t segStart = 0, segEnd = 0, prev = 0;
    bool open = false;

    auto closeSegment = [&]() {
        if (open && segEnd > segStart) segments.push_back({segStart, segEnd});
        open = false;
    };

    while (cursor < runs.size() && runs[cursor].start + runs[cursor].len <= from) ++cursor;
    for (size_t i = cursor; i < runs.size(); ++i) {
        const CoverRun& run = runs[i];
        uint64_t start = std::max(from, run.start);
        uint64_t end = std::min(to, run.start + run.len);
        if (run.start >= to) break;
        if (end <= start) continue;

        if (open) {
            if ((start - prev) > maxBlockDist) {
                closeSegment();
            } else {
                cumulative -= weight * static_cast<double>(start - prev);
                if (cumulative < 0.0) closeSegment();
            }
        }
        if (!open) {
            open = true;
            segStart = start;
            segEnd = start;
            cumulative = 0.0;
        }
        cumulative += static_cast<double>(end - start);
        if (cumulative >= 0.0) segEnd = end;
        prev = end;
    }
    closeSegment();
}

// junction class of each interstitial row against its nearest row, within link; rows are start-ascending
void assignJunctions(std::vector<TelomereBlock*>& rows, uint32_t link) {
    for (size_t i = 0; i < rows.size(); ++i) {
        if (rows[i]->pieces > 0) continue;
        uint64_t prevEnd = (i > 0) ? rows[i - 1]->start + rows[i - 1]->blockLen : 0;
        uint64_t nextStart = (i + 1 < rows.size()) ? rows[i + 1]->start : 0;
        uint64_t curEnd = rows[i]->start + rows[i]->blockLen;
        uint64_t distPrev = (i > 0) ? (rows[i]->start > prevEnd ? rows[i]->start - prevEnd : 0) : UINT64_MAX;
        uint64_t distNext = (i + 1 < rows.size())
            ? (nextStart > curEnd ? nextStart - curEnd : 0) : UINT64_MAX;
        bool useNext = distNext < distPrev;
        uint64_t dist = useNext ? distNext : distPrev;
        if (dist > link) {
            rows[i]->junction = 's';
            continue;
        }
        char first  = useNext ? rows[i]->strandLabel     : rows[i - 1]->strandLabel;
        char second = useNext ? rows[i + 1]->strandLabel : rows[i]->strandLabel;
        if (first == 'b' || second == 'b') rows[i]->junction = 's';
        else if (first == 'q' && second == 'p') rows[i]->junction = 'f';
        else if (first == 'p' && second == 'q') rows[i]->junction = 't';
        else rows[i]->junction = 'r';
    }
}
} // namespace


// per contig: chain qualifying same-strand pieces inward from the anchored end, stopping at a real opposite array
TelomereBlock Teloscope::getTerminalBlocks(const std::vector<MatchInfo>& matches,
        const std::vector<CoverRun>& runsFwd, const std::vector<CoverRun>& runsRev,
        uint64_t contigStart, uint64_t contigEnd, bool fromStart,
        std::vector<std::pair<uint64_t, uint64_t>>& outPieces, uint64_t& outProbe) {
    TelomereBlock chain;

    const uint32_t maxBlockDist = userInput.maxBlockDist;
    const uint32_t minLen = userInput.minBlockLen;
    const uint32_t minCounts = userInput.minBlockCounts;
    const uint32_t link = userInput.linkDistance;
    const float density = userInput.minBlockDensity;
    const float weight = (density >= 1.0f) ? 1e9f : density / (1.0f - density);
    const uint32_t zone = std::min(userInput.terminalTolerance, userInput.terminalLimit);

    // R6: furthest coordinate any trim looked toward the interior, plus its bridging margin
    uint64_t probeBound = fromStart ? contigStart : contigEnd;
    auto trackProbe = [&](uint64_t probe) {
        probeBound = fromStart ? std::max(probeBound, probe + maxBlockDist)
                                : std::min(probeBound, probe > maxBlockDist ? probe - maxBlockDist : 0);
    };

    // an array is real when its own unclamped trim reaches -l, -c and is strand-pure (R1)
    auto isRealArray = [&](const std::vector<CoverRun>& runs, uint64_t anchorEdge) -> bool {
        uint64_t farEdge = fromStart ? contigEnd : contigStart;
        uint64_t probe;
        uint64_t trimmed = trimInward(runs, anchorEdge, farEdge, maxBlockDist, weight, fromStart, probe);
        trackProbe(probe);
        uint64_t lo = fromStart ? anchorEdge : trimmed, hi = fromStart ? trimmed : anchorEdge;
        if (hi <= lo || (hi - lo) < minLen) return false;
        TelomereBlock tally;
        tally.start = lo;
        tally.blockLen = static_cast<uint32_t>(hi - lo);
        tallyBlock(matches, tally);
        bool isFwd = (&runs == &runsFwd);
        uint32_t count = isFwd ? tally.fwdCanCount : tally.revCanCount;
        char label = computeStrandLabel(tally.fwdCanCount, tally.fwdCanCount + tally.revCanCount, userInput.labelThreshold);
        return count >= minCounts && label == (isFwd ? 'p' : 'q');
    };

    // trim one piece from its anchor, stopping at the first real opposite-strand array
    auto trimPiece = [&](uint64_t anchorEdge, bool fwd, uint64_t farEdge) -> uint64_t {
        const std::vector<CoverRun>& runsSame = fwd ? runsFwd : runsRev;
        const std::vector<CoverRun>& runsOpp  = fwd ? runsRev : runsFwd;
        uint64_t probe;
        uint64_t trimmed = trimInward(runsSame, anchorEdge, farEdge, maxBlockDist, weight, fromStart, probe);
        trackProbe(probe);

        if (fromStart) {
            for (const CoverRun& r : runsOpp) {
                if (r.start <= anchorEdge) continue;
                if (r.start >= trimmed) break;
                if (isRealArray(runsOpp, r.start)) {
                    trimmed = trimInward(runsSame, anchorEdge, r.start, maxBlockDist, weight, fromStart, probe);
                    trackProbe(probe);
                    break;
                }
            }
        } else {
            for (auto it = runsOpp.rbegin(); it != runsOpp.rend(); ++it) {
                uint64_t rEnd = it->start + it->len;
                if (rEnd >= anchorEdge) continue;
                if (rEnd <= trimmed) break;
                if (isRealArray(runsOpp, rEnd)) {
                    trimmed = trimInward(runsSame, anchorEdge, rEnd, maxBlockDist, weight, fromStart, probe);
                    trackProbe(probe);
                    break;
                }
            }
        }
        return trimmed;
    };

    std::vector<MatchInfo> canon;
    canon.reserve(matches.size());
    for (const MatchInfo& m : matches) if (m.isCanonical) canon.push_back(m);
    size_t n = canon.size();
    if (n == 0) { outProbe = probeBound; return chain; }

    uint64_t reach = fromStart ? std::min(contigEnd, contigStart + zone)
                                : (contigEnd > zone ? std::max(contigStart, contigEnd - zone) : contigStart);
    bool haveStrand = false;
    bool curFwd = false;
    std::vector<std::pair<uint64_t, uint64_t>> pieces;

    for (size_t k = 0; k < n; ++k) {
        const MatchInfo& m = fromStart ? canon[k] : canon[n - 1 - k];
        uint64_t edge = fromStart ? static_cast<uint64_t>(m.position)
                                   : static_cast<uint64_t>(m.position) + m.matchSize;
        if (fromStart ? (edge > reach) : (edge < reach)) break; // a match exactly at reach is still within link

        bool mFwd = m.isForward;
        if (haveStrand && mFwd != curFwd) {
            if (isRealArray(mFwd ? runsFwd : runsRev, edge)) break;
            continue;
        }

        uint64_t farEdge = fromStart ? contigEnd : contigStart;
        uint64_t trimmed = trimPiece(edge, mFwd, farEdge);
        uint64_t pStart = fromStart ? edge : trimmed;
        uint64_t pEnd = fromStart ? trimmed : edge;
        if (pEnd <= pStart || (pEnd - pStart) < minLen) continue;

        TelomereBlock tally;
        tally.start = pStart;
        tally.blockLen = static_cast<uint32_t>(pEnd - pStart);
        tallyBlock(matches, tally);
        uint32_t canCount = mFwd ? tally.fwdCanCount : tally.revCanCount;
        if (canCount < minCounts) continue;
        char pieceLabel = computeStrandLabel(tally.fwdCanCount, tally.fwdCanCount + tally.revCanCount, userInput.labelThreshold);
        if (pieceLabel != (mFwd ? 'p' : 'q')) continue; // pieces are strand-pure (R1)

        pieces.push_back({pStart, pEnd});
        chain.fwdCanCount += tally.fwdCanCount;
        chain.revCanCount += tally.revCanCount;
        chain.fwdNonCanCount += tally.fwdNonCanCount;
        chain.revNonCanCount += tally.revNonCanCount;
        chain.teloLen += (pEnd - pStart);
        curFwd = mFwd;
        haveStrand = true;
        reach = fromStart ? (pEnd + link) : (pStart > link ? pStart - link : 0);

        // past the last piece: skip everything it already consumed
        while (k + 1 < n) {
            const MatchInfo& m2 = fromStart ? canon[k + 1] : canon[n - 1 - (k + 1)];
            uint64_t edge2 = fromStart ? static_cast<uint64_t>(m2.position)
                                        : static_cast<uint64_t>(m2.position) + m2.matchSize;
            bool inside = fromStart ? (edge2 < pEnd) : (edge2 > pStart);
            if (!inside) break;
            ++k;
        }
    }

    // R6: the furthest of every trim's reach and the final anchor-search bound
    outProbe = fromStart ? std::max(probeBound, reach) : std::min(probeBound, reach);
    if (pieces.empty()) return chain;
    uint64_t lo = pieces.front().first, hi = pieces.front().second;
    for (const auto& p : pieces) {
        lo = std::min(lo, p.first);
        hi = std::max(hi, p.second);
    }
    chain.start = lo;
    chain.blockLen = static_cast<uint32_t>(hi - lo);
    chain.pieces = static_cast<uint16_t>(pieces.size());
    chain.strandLabel = curFwd ? 'p' : 'q';
    chain.anchorSide = fromStart ? 'p' : 'q';
    outPieces = std::move(pieces);
    return chain;
}
// outside the terminal chains' spans: -k seeds with the inversion cut, grouped within -d, trimmed by getMaximalSegments
void Teloscope::getInterstitialBlocks(const std::vector<MatchInfo>& matches,
        const std::vector<CoverRun>& allRuns, uint64_t regionStart, uint64_t regionEnd,
        std::vector<TelomereBlock>& outBlocks) {
    if (regionEnd <= regionStart) return;

    auto lo = std::lower_bound(matches.begin(), matches.end(), regionStart,
        [](const MatchInfo& m, uint64_t v) { return m.position < v; });
    auto hi = std::lower_bound(matches.begin(), matches.end(), regionEnd,
        [](const MatchInfo& m, uint64_t v) { return m.position < v; });

    std::vector<Seed> seeds;
    getSeeds(lo, hi, userInput.maxMatchDist, seeds);
    if (seeds.empty()) return;

    const float density = userInput.minBlockDensity;
    const float weight = (density >= 1.0f) ? 1e9f : density / (1.0f - density);
    const uint32_t maxBlockDist = userInput.maxBlockDist;
    size_t segCursor = 0, tallyCursor = 0, coverCursor = 0; // groups/segments are start-ascending

    size_t i = 0;
    while (i < seeds.size()) {
        size_t j = i;
        uint64_t groupStart = seeds[i].start, groupEnd = seeds[i].end;
        bool fwd = seeds[i].isForward;
        while (j + 1 < seeds.size() && seeds[j + 1].isForward == fwd &&
               seeds[j + 1].start <= groupEnd + maxBlockDist) {
            ++j;
            groupEnd = std::max(groupEnd, seeds[j].end);
        }
        groupStart = std::max(groupStart, regionStart);
        groupEnd = std::min(groupEnd, regionEnd);

        std::vector<std::pair<uint64_t, uint64_t>> segments;
        if (groupEnd > groupStart) {
            getMaximalSegments(allRuns, groupStart, groupEnd, maxBlockDist, weight, segCursor, segments);
        }

        for (const auto& seg : segments) {
            TelomereBlock block;
            block.start = seg.first;
            block.blockLen = static_cast<uint32_t>(seg.second - seg.first);
            tallyBlock(matches, block, tallyCursor);
            const uint64_t canonicalCount =
                static_cast<uint64_t>(block.fwdCanCount) + block.revCanCount;
            if (canonicalCount < minCanonicalCount) continue;
            if (getCoveredBases(allRuns, seg.first, seg.second, coverCursor) < density * block.blockLen) continue;

            const uint64_t forwardCount =
                static_cast<uint64_t>(block.fwdCanCount) + block.fwdNonCanCount;
            const uint64_t totalCount = forwardCount + block.revCanCount + block.revNonCanCount;
            block.strandLabel = computeStrandLabel(forwardCount, totalCount, userInput.labelThreshold);
            block.teloLen = block.blockLen;
            outBlocks.push_back(block);
        }
        i = j + 1;
    }
}


void Teloscope::labelTerminalBlocks(
    std::vector<TelomereBlock>& blocks,
    std::string& terminalLabel, ScaffoldType& scaffoldType, uint8_t& anomalyFlags) {

    terminalLabel.clear();
    anomalyFlags = 0;

    std::sort(blocks.begin(), blocks.end(),
            [](const TelomereBlock &a, const TelomereBlock &b) {
                return a.start < b.start;
            });

    TelomereBlock* armP = nullptr;
    TelomereBlock* armQ = nullptr;

    // uppercase for arms, lowercase for contig rows, '*' after a discordant row
    for (auto& block : blocks) {
        char letter = block.isScaffold ? static_cast<char>(std::toupper(block.anchorSide))
                                        : block.anchorSide;
        terminalLabel += letter;
        if (block.strandLabel != block.anchorSide) terminalLabel += '*';

        if (!block.isScaffold) continue;
        if (block.anchorSide == 'p') armP = &block; else armQ = &block;
    }

    bool hasP = (armP != nullptr);
    bool hasQ = (armQ != nullptr);

    if (hasP && hasQ) scaffoldType = ScaffoldType::T2T;
    else if (hasP || hasQ) scaffoldType = ScaffoldType::INCOMPLETE;
    else scaffoldType = ScaffoldType::NONE;

    if (hasP) {
        if (armP->strandLabel != armP->anchorSide) anomalyFlags |= ANOMALY_DISCORDANT_P;
        if (armP->pieces >= 2) anomalyFlags |= ANOMALY_FRAGMENTED_P;
    }
    if (hasQ) {
        if (armQ->strandLabel != armQ->anchorSide) anomalyFlags |= ANOMALY_DISCORDANT_Q;
        if (armQ->pieces >= 2) anomalyFlags |= ANOMALY_FRAGMENTED_Q;
    }
}


void Teloscope::analyzeWindow(const std::string_view &window, uint64_t windowStart,
                            WindowData& windowData, WindowData& nextOverlapData,
                            SegmentData& segmentData, uint64_t segmentSize, uint64_t absPos) {

    windowData.windowStart = windowStart;
    unsigned short int longestPatternSize = this->trie.getLongestPatternSize();
    uint32_t step = userInput.step;
    uint32_t overlapSize = userInput.windowSize - step;
    uint32_t terminalLimit = userInput.terminalLimit;
    bool computeGC = userInput.outGC;
    bool computeEntropy = userInput.outEntropy;
    bool hasLastCanonical = false;
    uint64_t lastCanonicalPos = 0;
    bool needMatchSeq = userInput.outMatches;

    // terminal status
    uint64_t windowEnd = windowStart + window.size();
    uint64_t terminalEnd = (segmentSize > terminalLimit) ? (segmentSize - terminalLimit) : 0;
    bool windowFullyTerminal = (windowEnd <= terminalLimit) || (windowStart >= terminalEnd);
    bool windowFullyInterstitial = (windowStart > terminalLimit) && (windowEnd < terminalEnd);

    // overlap conditions
    bool alwaysMainWindow = (overlapSize == 0 || windowStart == 0);
    bool hasOverlap = (overlapSize != 0);

    // pushes are partitioned by match start; the scan keeps its lookback so a
    // canonical dimer straddling the overlap edge is still seen
    uint32_t startIndex = (alwaysMainWindow || std::min(step, overlapSize) <= longestPatternSize)
                            ? 0
                            : std::min(step, overlapSize) - longestPatternSize;

    for (uint32_t i = startIndex; i < window.size(); ++i) {
        if (computeGC || computeEntropy) {
            char nucleotide = window[i];
            uint8_t index;
            switch (nucleotide) {
                case 'A': index = 0; break;
                case 'C': index = 1; break;
                case 'G': index = 2; break;
                case 'T': index = 3; break;
                default: continue;
            }

            if (alwaysMainWindow || i >= overlapSize) {
                windowData.nucleotideCounts[index]++;
            }
            if (hasOverlap && i >= step) {
                nextOverlapData.nucleotideCounts[index]++;
            }
        }

        // trie pattern matching
        int32_t current = trie.getRoot();
        // probe past the window end so a match straddling a boundary is not lost
        uint32_t scanLimit = static_cast<uint32_t>(
            std::min<uint64_t>(static_cast<uint64_t>(i) + longestPatternSize, segmentSize - windowStart));

        for (uint32_t j = i; j < scanLimit; ++j) {
            current = trie.getChild(current, window.data()[j]);
            if (current < 0) break;

            if (trie.isEnd(current)) {
                uint16_t matchLen = static_cast<uint16_t>(j - i + 1);
                bool isForward = trie.isForward(current);
                bool isCanonical = trie.isCanonical(current);
                uint64_t matchPos = absPos + windowStart + i;

                bool isTerminal;
                if (windowFullyTerminal) {
                    isTerminal = true;
                } else if (windowFullyInterstitial) {
                    isTerminal = false;
                } else {
                    uint64_t absI = windowStart + i;
                    isTerminal = (absI <= terminalLimit || absI >= terminalEnd);
                }

                MatchInfo matchInfo{matchPos, matchLen, isCanonical, isForward};

                // Check dimers
                if (isCanonical) {
                    if (hasLastCanonical && (matchPos - lastCanonicalPos) <= userInput.canonicalSize) {
                        if (!windowData.hasCanDimer && (alwaysMainWindow || i >= overlapSize)) {
                            windowData.hasCanDimer = true;
                        }
                        if (!nextOverlapData.hasCanDimer && hasOverlap && i >= step) {
                            nextOverlapData.hasCanDimer = true;
                        }
                    }
                    lastCanonicalPos = matchPos;
                    hasLastCanonical = true;
                }

                // Update windowData
                if (alwaysMainWindow || i >= overlapSize) {
                    if (isCanonical) {
                        windowData.canonicalCounts++;
                        windowData.canonicalCovered += matchLen;
                        segmentData.canonicalCounts++;
                        if (needMatchSeq) {
                            segmentData.canonicalMatches.push_back({matchPos, std::string(window.data() + i, matchLen)});
                        }
                    } else {
                        windowData.nonCanonicalCounts++;
                        windowData.nonCanonicalCovered += matchLen;
                        if (isTerminal && needMatchSeq) {
                            segmentData.nonCanonicalMatches.push_back({matchPos, std::string(window.data() + i, matchLen)});
                        }
                    }

                    // route by strand
                    if (isForward) {
                        windowData.fwdCounts++;
                        windowData.fwdCovered += matchLen;
                        segmentData.fwdCounts++;
                    } else {
                        windowData.revCounts++;
                        windowData.revCovered += matchLen;
                        segmentData.revCounts++;
                    }
                    segmentData.allMatches.push_back(matchInfo);
                }

                // Update nextOverlapData
                if (hasOverlap && i >= step) {
                    if (isCanonical) {
                        nextOverlapData.canonicalCounts++;
                        nextOverlapData.canonicalCovered += matchLen;
                    } else {
                        nextOverlapData.nonCanonicalCounts++;
                        nextOverlapData.nonCanonicalCovered += matchLen;
                    }

                    // fwd/rev overlap metrics
                    if (isForward) {
                        nextOverlapData.fwdCounts++;
                        nextOverlapData.fwdCovered += matchLen;
                    } else {
                        nextOverlapData.revCounts++;
                        nextOverlapData.revCovered += matchLen;
                    }
                }

            }
        }
    }
}


SegmentData Teloscope::scanSegment(std::string &sequence, uint64_t absPos,
                                   bool isFirst, bool isLast, bool buildP, bool buildQ) {
    SegmentData segmentData;
    uint64_t segmentSize = sequence.size();
    if (absPos + segmentSize > (1ULL << 40)) {
        // a worker thread cannot exit through the pool it is still parked in
        std::cerr << "Error: sequence coordinate exceeds the 1.1 Tb limit." << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    unsigned short int longestPatternSize = this->trie.getLongestPatternSize();
    bool useWindowed = userInput.outWinRepeats || userInput.outGC ||
                       userInput.outEntropy || userInput.outMatches;

    // fast mode only: the tiled head/tail window edges, for the interstitial margin filter below
    bool fastTiled = false, fastHeadActive = false, fastTailActive = false;
    uint64_t fastHeadEnd = 0, fastTailStart = 0;

    if (!useWindowed) {
        // ========== cheap path: trie matching only, no window bookkeeping ==========
        auto processRegion = [&](uint64_t start, uint64_t end, std::vector<MatchInfo>& out) {
            for (uint64_t i = start; i < end; ++i) {
                int32_t node = trie.getRoot();
                uint64_t scanLimit = std::min(i + static_cast<uint64_t>(longestPatternSize), segmentSize);

                for (uint64_t j = i; j < scanLimit; ++j) {
                    node = trie.getChild(node, sequence[j]);
                    if (node < 0) break;

                    if (trie.isEnd(node)) {
                        uint16_t len = static_cast<uint16_t>(j - i + 1);
                        bool isForward = trie.isForward(node);
                        bool isCanonical = trie.isCanonical(node);

                        out.push_back({absPos + i, len, isCanonical, isForward});
                        if (isCanonical) segmentData.canonicalCounts++;
                        if (isForward) segmentData.fwdCounts++;
                        else segmentData.revCounts++;
                    }
                }
            }
        };

        uint32_t t = userInput.terminalLimit;
        bool wholeContig = !userInput.ultraFastMode ||
                           (buildP && buildQ && segmentSize <= 2ULL * t);

        if (wholeContig) {
            processRegion(0, segmentSize, segmentData.allMatches);
        } else {
            // R6: tiled, non-overlapping ranges extended while the terminal builder still looks that far
            uint64_t h  = buildP ? std::min<uint64_t>(t, segmentSize) : 0;
            uint64_t t0 = buildQ ? (segmentSize > t ? segmentSize - t : 0) : segmentSize;
            std::vector<MatchInfo> headMatches, tailMatches;
            if (buildP) processRegion(0, h, headMatches);
            if (buildQ) processRegion(t0, segmentSize, tailMatches);

            bool grew = true;
            while (grew) {
                grew = false;
                if (buildP && h < t0) {
                    std::vector<CoverRun> hFwd, hRev;
                    getCoverRuns(headMatches, Orient::Fwd, hFwd);
                    getCoverRuns(headMatches, Orient::Rev, hRev);
                    std::vector<std::pair<uint64_t, uint64_t>> tmpPieces;
                    uint64_t probeP;
                    getTerminalBlocks(headMatches, hFwd, hRev, absPos, absPos + segmentSize, true, tmpPieces, probeP);
                    if (probeP + longestPatternSize >= absPos + h) {
                        uint64_t newH = std::min<uint64_t>(h + t, t0);
                        processRegion(h, newH, headMatches);
                        h = newH;
                        grew = true;
                    }
                }
                if (buildQ && t0 > h) {
                    std::vector<CoverRun> tFwd, tRev;
                    getCoverRuns(tailMatches, Orient::Fwd, tFwd);
                    getCoverRuns(tailMatches, Orient::Rev, tRev);
                    std::vector<std::pair<uint64_t, uint64_t>> tmpPieces;
                    uint64_t probeQ;
                    getTerminalBlocks(tailMatches, tFwd, tRev, absPos, absPos + segmentSize, false, tmpPieces, probeQ);
                    if (probeQ <= absPos + t0 + longestPatternSize) {
                        uint64_t newT0 = (t0 > t) ? std::max<uint64_t>(t0 - t, h) : h;
                        std::vector<MatchInfo> extra;
                        processRegion(newT0, t0, extra);
                        extra.insert(extra.end(), tailMatches.begin(), tailMatches.end());
                        tailMatches.swap(extra);
                        t0 = newT0;
                        grew = true;
                    }
                }
            }
            fastTiled = true;
            fastHeadActive = buildP;
            fastTailActive = buildQ;
            fastHeadEnd = absPos + h;
            fastTailStart = absPos + t0;
            segmentData.allMatches = std::move(headMatches);
            segmentData.allMatches.insert(segmentData.allMatches.end(), tailMatches.begin(), tailMatches.end());
        }

    } else {
        // ========== full path: window-based scan (always whole-contig) ==========
        uint32_t windowSize = userInput.windowSize;
        uint32_t step = userInput.step;

        // capped reserve
        constexpr uint64_t maxMatchReserve = 1000000;
        if (userInput.outMatches) {
            segmentData.canonicalMatches.reserve(std::min(segmentSize / 6, maxMatchReserve));
            segmentData.nonCanonicalMatches.reserve(std::min(segmentSize / 6, maxMatchReserve));
        }
        segmentData.allMatches.reserve(std::min(segmentSize / 3, maxMatchReserve));

        bool keepWindows = userInput.outWinRepeats || userInput.outEntropy || userInput.outGC;
        if (keepWindows && segmentSize > windowSize) {
            segmentData.windows.reserve((segmentSize - windowSize) / step + 2);
        }

        WindowData prevOverlapData; // Data from previous overlap
        WindowData nextOverlapData; // Data for next overlap

        std::vector<WindowData> windows;
        uint64_t windowStart = 0;
        uint64_t currentWindowSize = std::min(static_cast<uint64_t>(windowSize), segmentSize);
        std::string_view windowView(sequence.data(), currentWindowSize);

        while (windowStart < segmentSize) {
            // Prepare and analyze current window
            WindowData windowData = prevOverlapData;

            analyzeWindow(windowView, windowStart,
                        windowData, nextOverlapData,
                        segmentData, segmentSize, absPos);

            if (userInput.outGC) {
                windowData.gcContent = getGCContent(windowData.nucleotideCounts, windowView.size());
            }
            if (userInput.outEntropy) {
                windowData.shannonEntropy = getShannonEntropy(windowData.nucleotideCounts, windowView.size());
            }

            // Update windowData
            windowData.windowStart = windowStart + absPos;
            windowData.currentWindowSize = currentWindowSize;
            segmentData.windowCounts++;
            if (keepWindows) windows.emplace_back(windowData);

            prevOverlapData = nextOverlapData;
            nextOverlapData = WindowData();

            windowStart += step;
            if (windowStart >= segmentSize) break;

            // Prepare next window
            currentWindowSize = std::min(static_cast<uint64_t>(windowSize), segmentSize - windowStart);
            windowView = std::string_view(sequence.data() + windowStart, currentWindowSize);
        }

        segmentData.windows = std::move(windows);
    }

    // ========== block building, per contig ==========
    std::vector<CoverRun> runsFwd, runsRev, allRuns;
    getCoverRuns(segmentData.allMatches, Orient::Fwd, runsFwd);
    getCoverRuns(segmentData.allMatches, Orient::Rev, runsRev);
    getCoverRuns(segmentData.allMatches, Orient::All, allRuns);

    uint64_t contigStart = absPos, contigEnd = absPos + segmentSize;
    TelomereBlock chainP, chainQ;
    std::vector<std::pair<uint64_t, uint64_t>> piecesP, piecesQ;
    uint64_t probeP, probeQ;
    if (buildP) chainP = getTerminalBlocks(segmentData.allMatches, runsFwd, runsRev, contigStart, contigEnd, true, piecesP, probeP);
    if (buildQ) chainQ = getTerminalBlocks(segmentData.allMatches, runsFwd, runsRev, contigStart, contigEnd, false, piecesQ, probeQ);

    bool haveP = chainP.pieces > 0, haveQ = chainQ.pieces > 0;
    if (haveP && haveQ && chainP.strandLabel == chainQ.strandLabel &&
        chainP.start + chainP.blockLen > chainQ.start) {
        // R3: same-strand overlap; merge the piece lists, not just the spans
        std::vector<std::pair<uint64_t, uint64_t>> merged = piecesP;
        merged.insert(merged.end(), piecesQ.begin(), piecesQ.end());
        std::sort(merged.begin(), merged.end());
        std::vector<std::pair<uint64_t, uint64_t>> mergedPieces;
        for (const auto& p : merged) {
            if (!mergedPieces.empty() && p.first <= mergedPieces.back().second) {
                mergedPieces.back().second = std::max(mergedPieces.back().second, p.second);
            } else {
                mergedPieces.push_back(p);
            }
        }
        TelomereBlock result;
        result.start = mergedPieces.front().first;
        result.blockLen = static_cast<uint32_t>(mergedPieces.back().second - result.start);
        result.pieces = static_cast<uint16_t>(mergedPieces.size());
        for (const auto& piece : mergedPieces) {
            TelomereBlock tally;
            tally.start = piece.first;
            tally.blockLen = static_cast<uint32_t>(piece.second - piece.first);
            tallyBlock(segmentData.allMatches, tally);
            result.teloLen += tally.blockLen;
            result.fwdCanCount += tally.fwdCanCount;
            result.revCanCount += tally.revCanCount;
            result.fwdNonCanCount += tally.fwdNonCanCount;
            result.revNonCanCount += tally.revNonCanCount;
        }
        result.strandLabel = chainP.strandLabel;
        result.anchorSide = result.strandLabel;
        result.isScaffold = (result.anchorSide == 'p') ? isFirst : isLast;
        chainP = result;
        haveQ = false;
    } else {
        if (haveP) chainP.isScaffold = isFirst;
        if (haveQ) chainQ.isScaffold = isLast;
    }

    uint64_t pFrom = haveP ? chainP.start : contigStart;
    uint64_t pTo   = haveP ? chainP.start + chainP.blockLen : contigStart;
    uint64_t qFrom = haveQ ? chainQ.start : contigEnd;
    uint64_t qTo   = haveQ ? chainQ.start + chainQ.blockLen : contigEnd;

    getInterstitialBlocks(segmentData.allMatches, allRuns, contigStart, pFrom, segmentData.interstitialBlocks);
    getInterstitialBlocks(segmentData.allMatches, allRuns, pTo, qFrom, segmentData.interstitialBlocks);
    getInterstitialBlocks(segmentData.allMatches, allRuns, qTo, contigEnd, segmentData.interstitialBlocks);

    if (fastTiled) {
        // fast mode: drop interstitial rows that a wider window could still extend, only at a real cut
        const uint64_t margin = userInput.maxBlockDist + userInput.maxMatchDist + longestPatternSize;
        bool headIsCut = fastHeadActive && fastHeadEnd != contigEnd;
        bool tailIsCut = fastTailActive && fastTailStart != contigStart;
        auto& rows = segmentData.interstitialBlocks;
        rows.erase(std::remove_if(rows.begin(), rows.end(), [&](const TelomereBlock& b) {
            bool nearHead = headIsCut && b.start < fastHeadEnd && (b.start + b.blockLen + margin >= fastHeadEnd);
            bool nearTail = tailIsCut && b.start >= fastTailStart && (b.start <= fastTailStart + margin);
            return nearHead || nearTail;
        }), rows.end());
    }

    std::vector<TelomereBlock*> contigRows;
    if (haveP) contigRows.push_back(&chainP);
    for (auto& block : segmentData.interstitialBlocks) contigRows.push_back(&block);
    if (haveQ) contigRows.push_back(&chainQ);
    std::sort(contigRows.begin(), contigRows.end(),
             [](const TelomereBlock* a, const TelomereBlock* b) { return a->start < b->start; });
    assignJunctions(contigRows, userInput.linkDistance);

    if (haveP) segmentData.terminalBlocks.push_back(chainP);
    if (haveQ) segmentData.terminalBlocks.push_back(chainQ);

    return segmentData;
}


void Teloscope::writeBEDFile(std::ofstream& windowDensityFile,
                            std::ofstream& windowCanonicalRatioFile,
                            std::ofstream& windowStrandRatioFile,
                            std::ofstream& windowGCFile,
                            std::ofstream& windowEntropyFile,
                            std::ofstream& canonicalMatchFile,
                            std::ofstream& noncanonicalMatchFile,
                            std::ofstream& terminalBlocksFile,
                            std::ofstream& interstitialBlocksFile,
                            std::ofstream& gapFile,
                            std::ofstream& reportFile) {

    // Keep BED and BEDGraph streams free of non-data comments. The always-written
    // report is the single run-level provenance record for their shared basename.
    writeProvenanceHeader(
        reportFile, userInput,
        userInput.ultraFastMode
            ? "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly\tgranular"
            : "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly\tgranular\tits\tcanonical\twindows");

    // BEDgraph headers
    if (userInput.outWinRepeats) {
        windowDensityFile << "track type=bedGraph name=\"Repeat Density\" description=\"Total repeat density per window\"\n";
        windowCanonicalRatioFile << "track type=bedGraph name=\"Canonical Ratio\" description=\"Canonical fraction of repeat density per window\"\n";
        windowStrandRatioFile << "track type=bedGraph name=\"Strand Ratio\" description=\"Forward-strand fraction of repeat density per window\"\n";
    }
    if (userInput.outEntropy) {
        windowEntropyFile << "track type=bedGraph name=\"Shannon Entropy\" description=\"Shannon entropy per window\"\n";
    }
    if (userInput.outGC) {
        windowGCFile << "track type=bedGraph name=\"GC Content\" description=\"GC content per window\"\n";
    }

    // Report header (console + file)
    std::cout << "\n+++ Path Summary Report +++\n";
    if (!userInput.ultraFastMode) {
        std::cout << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly\tgranular\tits\tcanonical\twindows\n";
        reportFile << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly\tgranular\tits\tcanonical\twindows\n";
    } else {
        std::cout << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly\tgranular\n";
        reportFile << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly\tgranular\n";
    }

    // Processing paths
    totalPaths = allPathData.size();
    std::vector<float> telomereLengths; // for getStats

    for (const auto& pathData : allPathData) {
        const auto& header = pathData.header;
        const auto& windows = pathData.windows;
        const auto& pos = pathData.seqPos;
        const uint16_t gaps = static_cast<uint16_t>(pathData.gapInfos.size());
        const auto& pathSize = pathData.pathSize;

        // Longest telomere blocks
        int longestCount = 0;
        std::string longestLabels;

        for (const auto& block : pathData.terminalBlocks) {
            if (block.isScaffold || userInput.manualCuration) {
                writeBlockRow(terminalBlocksFile, header, block, pathSize,
                             block.isScaffold ? "scaffold" : "contig");
            }

            if (block.isScaffold) {
                longestCount++;
                longestLabels += block.anchorSide;
                telomereLengths.push_back(static_cast<float>(block.teloLen));
            }
        }

        for (const auto& block : pathData.interstitialBlocks) {
            writeBlockRow(interstitialBlocksFile, header, block, pathSize,
                         junctionToString(block.junction));
        }

        // Gap positions
        for (const auto& gap : pathData.gapInfos) {
            gapFile << header << "\t"
                    << gap.start << "\t"
                    << (gap.start + gap.length) << "\n";
        }

        // All canonical and terminal non-canonical matches
        if (userInput.outMatches) {
            for (const auto& match : pathData.canonicalMatches) {
                canonicalMatchFile << header << "\t"
                                << match.position << "\t"
                                << (match.position + match.matchSeq.size()) << "\t"
                                << match.matchSeq << "\n";
            }

            for (const auto& match : pathData.nonCanonicalMatches) {
                noncanonicalMatchFile << header << "\t"
                                    << match.position << "\t"
                                    << (match.position + match.matchSeq.size()) << "\t"
                                    << match.matchSeq << "\n";
            }
        }

        // Process window metrics
        for (const auto& window : windows) {
            uint64_t windowEnd = window.windowStart + window.currentWindowSize;

            if (userInput.outWinRepeats) {
                uint32_t totalCovered = window.fwdCovered + window.revCovered;
                float totalDensity = static_cast<float>(totalCovered) / window.currentWindowSize;
                float canonRatio = (totalCovered > 0)
                    ? static_cast<float>(window.canonicalCovered) / (window.canonicalCovered + window.nonCanonicalCovered)
                    : -1.0f;
                float strandRatio = (totalCovered > 0)
                    ? static_cast<float>(window.fwdCovered) / (window.fwdCovered + window.revCovered)
                    : -1.0f;

                windowDensityFile << header << "\t" << window.windowStart << "\t" << windowEnd
                                    << "\t" << totalDensity << "\n";
                windowCanonicalRatioFile << header << "\t" << window.windowStart << "\t" << windowEnd
                                           << "\t" << canonRatio << "\n";
                windowStrandRatioFile << header << "\t" << window.windowStart << "\t" << windowEnd
                                        << "\t" << strandRatio << "\n";
            }
            if (userInput.outEntropy) {
                windowEntropyFile << header << "\t" << window.windowStart << "\t" << windowEnd
                                    << "\t" << window.shannonEntropy << "\n";
            }
            if (userInput.outGC) {
                windowGCFile << header << "\t" << window.windowStart << "\t" << windowEnd
                               << "\t" << window.gcContent << "\n";
            }
        }

        // Output path summary (console + file)
        const char* typeStr = scaffoldTypeToString(pathData.scaffoldType);
        const std::string anomalyStr = anomalyFlagsToString(pathData.anomalyFlags);
        const char* labelsStr = longestLabels.empty() ? "none" : longestLabels.c_str();

        std::cout << pos + 1 << "\t" << header << "\t"
                << longestCount << "\t" << labelsStr << "\t"
                << gaps << "\t" << typeStr << "\t"
                << anomalyStr << "\t"
                << pathData.terminalLabel;
        reportFile << pos + 1 << "\t" << header << "\t"
                << longestCount << "\t" << labelsStr << "\t"
                << gaps << "\t" << typeStr << "\t"
                << anomalyStr << "\t"
                << pathData.terminalLabel;

        totalTelomeres += longestCount;
        if (longestCount == 0) totalZeroTelomeres++;
        else if (longestCount == 1) totalOneTelomere++;
        else totalTwoTelomeres++;
        totalGaps += gaps;

        // Expand path summary
        if (!userInput.ultraFastMode) {
            std::cout << "\t"
                    << pathData.interstitialBlocks.size() << "\t"
                    << pathData.canonicalCounts << "\t"
                    << pathData.windowCounts;
            reportFile << "\t"
                    << pathData.interstitialBlocks.size() << "\t"
                    << pathData.canonicalCounts << "\t"
                    << pathData.windowCounts;

            totalNWindows += pathData.windowCounts;
            totalITS += pathData.interstitialBlocks.size();
            totalCanMatches += pathData.canonicalCounts;
        }
        std::cout << "\n";
        reportFile << "\n";
    }

    // Calculate telomere statistics
    if (!telomereLengths.empty()) {
        Stats stats = getStats(telomereLengths);
        teloMean = stats.mean;
        teloMedian = stats.median;
        teloMin = stats.min;
        teloMax = stats.max;
    }

}


void Teloscope::handleBEDFile() {
    lg.verbose("\nReporting window matches and metrics...");

    constexpr size_t ioBufSize = 1 << 20; // 1MB write buffer per file

    // buffers must outlive the ofstreams that reference them
    std::vector<char> densityBuf(ioBufSize), canonRatioBuf(ioBufSize), strandRatioBuf(ioBufSize);
    std::vector<char> gcBuf(ioBufSize), entropyBuf(ioBufSize);
    std::vector<char> canonMatchBuf(ioBufSize), noncanonMatchBuf(ioBufSize);
    std::vector<char> termBlockBuf(ioBufSize), itsBlockBuf(ioBufSize), gapBuf(ioBufSize);
    std::vector<char> reportBuf(ioBufSize);

    std::ofstream windowDensityFile;
    std::ofstream windowCanonicalRatioFile;
    std::ofstream windowStrandRatioFile;
    std::ofstream windowGCFile;
    std::ofstream windowEntropyFile;
    std::ofstream canonicalMatchFile;
    std::ofstream noncanonicalMatchFile;
    std::ofstream terminalBlocksFile;
    std::ofstream interstitialBlocksFile;
    std::ofstream gapFile;
    std::ofstream reportFile;

    std::string base = userInput.outRoute + "/" + userInput.inSequenceName;

    auto openFile = [&](std::ofstream& file, const std::string& path, std::vector<char>& buf) {
        file.rdbuf()->pubsetbuf(buf.data(), ioBufSize);
        file.open(path);
        if (!file.is_open()) {
            fprintf(stderr, "Error: Could not open '%s' for writing.\n", path.c_str());
            exit(EXIT_FAILURE);
        }
    };

    if (userInput.outWinRepeats) {
        openFile(windowDensityFile, base + "_window_repeat_density.bedgraph", densityBuf);
        openFile(windowCanonicalRatioFile, base + "_window_canonical_ratio.bedgraph", canonRatioBuf);
        openFile(windowStrandRatioFile, base + "_window_strand_ratio.bedgraph", strandRatioBuf);
    }
    if (userInput.outGC) {
        openFile(windowGCFile, base + "_window_gc.bedgraph", gcBuf);
    }
    if (userInput.outEntropy) {
        openFile(windowEntropyFile, base + "_window_entropy.bedgraph", entropyBuf);
    }

    if (userInput.outMatches) {
        openFile(canonicalMatchFile, base + "_canonical_matches.bed", canonMatchBuf);
        openFile(noncanonicalMatchFile, base + "_noncanonical_matches.bed", noncanonMatchBuf);
    }

    openFile(interstitialBlocksFile, base + "_interstitial_telomeres.bed", itsBlockBuf);
    openFile(terminalBlocksFile, base + "_terminal_telomeres.bed", termBlockBuf);
    openFile(gapFile, base + "_gaps.bed", gapBuf);
    openFile(reportFile, base + "_report.tsv", reportBuf);

    writeBEDFile(windowDensityFile, windowCanonicalRatioFile,
                windowStrandRatioFile,
                windowGCFile, windowEntropyFile,
                canonicalMatchFile, noncanonicalMatchFile,
                terminalBlocksFile, interstitialBlocksFile,
                gapFile, reportFile);

    printSummary(reportFile);
    reportFile.close();

    // Close all files
    if (userInput.outWinRepeats) {
        windowDensityFile.close();
        windowCanonicalRatioFile.close();
        windowStrandRatioFile.close();
    }
    if (userInput.outGC) {
        windowGCFile.close();
    }
    if (userInput.outEntropy) {
        windowEntropyFile.close();
    }

    if (userInput.outMatches) {
        canonicalMatchFile.close();
        noncanonicalMatchFile.close();
    }

    interstitialBlocksFile.close();
    terminalBlocksFile.close();
    gapFile.close();
}


void Teloscope::computeSummaryCounts() {
    std::vector<uint64_t> scaffoldLens;
    std::vector<uint64_t> contigLens;
    scaffoldLens.reserve(allPathData.size());

    for (const auto& pathData : allPathData) {
        const bool hasGaps = !pathData.gapInfos.empty();
        switch (pathData.scaffoldType) {
            case ScaffoldType::T2T:        hasGaps ? totalGappedT2T++ : totalT2T++; break;
            case ScaffoldType::INCOMPLETE: hasGaps ? totalGappedIncomplete++ : totalIncomplete++; break;
            case ScaffoldType::NONE:       hasGaps ? totalGappedNone++ : totalNone++; break;
        }

        const uint8_t flags = pathData.anomalyFlags;
        if (flags) totalFlagged++;
        if (flags & ANOMALY_DISCORDANT_P) totalDiscordantArms++;
        if (flags & ANOMALY_DISCORDANT_Q) totalDiscordantArms++;
        if (flags & ANOMALY_FRAGMENTED_P) totalFragmented++;
        if (flags & ANOMALY_FRAGMENTED_Q) totalFragmented++;

        // contig lengths = runs between gaps
        scaffoldLens.push_back(pathData.pathSize);
        std::vector<GapInfo> gaps = pathData.gapInfos;
        std::sort(gaps.begin(), gaps.end(),
                  [](const GapInfo& a, const GapInfo& b) { return a.start < b.start; });
        uint64_t prevEnd = 0;
        for (const auto& g : gaps) {
            if (g.start > prevEnd) contigLens.push_back(g.start - prevEnd);
            prevEnd = g.start + g.length;
        }
        if (pathData.pathSize > prevEnd) contigLens.push_back(pathData.pathSize - prevEnd);
    }

    scaffoldN50 = computeN50(scaffoldLens);
    contigN50 = computeN50(contigLens);
}


void Teloscope::printSummary(std::ofstream& reportFile) {
    computeSummaryCounts();

    auto out = [&](const auto&... args) {
        std::ostringstream ss;
        (ss << ... << args);
        std::string s = ss.str();
        std::cout << s;
        reportFile << s;
    };

    out("\n+++ Assembly Summary Report +++\n");
    out("Total paths:\t", totalPaths, "\n");
    if (userInput.sequenceFilterActive) {
        out("Filter input paths:\t", userInput.filterInputCount, "\n");
        out("Filter selected paths:\t", userInput.filterSelectedCount, "\n");
    }
    out("Total gaps:\t", totalGaps, "\n");
    out("Scaffold N50:\t", scaffoldN50, "\n");
    out("Contig N50:\t", contigN50, "\n");
    out("Total telomeres:\t", totalTelomeres, "\n");

    if (!userInput.ultraFastMode) {
        out("Total ITS blocks:\t", totalITS, "\n");
        out("Total canonical matches:\t", totalCanMatches, "\n");
        out("Total windows analyzed:\t", totalNWindows, "\n");
    }

    out("\n+++ Telomere Statistics +++\n");
    if (totalTelomeres > 0) {
        out("Mean length:\t", teloMean, "\n");
        out("Median length:\t", teloMedian, "\n");
        out("Min length:\t", teloMin, "\n");
        out("Max length:\t", teloMax, "\n");
    }
    else {
        out("No telomeres found for statistics.\n");
    }

    out("\n+++ Chromosome Telomere Counts+++\n");
    out("Two telomeres:\t", totalTwoTelomeres, "\n");
    out("One telomere:\t", totalOneTelomere, "\n");
    out("Zero telomeres:\t", totalZeroTelomeres, "\n");

    // these six partition the scaffolds, which is the property worth protecting
    out("\n+++ Chromosome Telomere/Gap Completeness+++\n");
    out("T2T:\t", totalT2T, "\n");
    out("Gapped T2T:\t", totalGappedT2T, "\n");

    out("Incomplete:\t", totalIncomplete, "\n");
    out("Gapped incomplete:\t", totalGappedIncomplete, "\n");

    out("No telomeres:\t", totalNone, "\n");
    out("Gapped no telomeres:\t", totalGappedNone, "\n");

    // plausibility, reported beside completeness. The detail lines may sum above
    // the flagged count, because one scaffold can carry more than one anomaly.
    out("\n+++ Scaffold Anomalies +++\n");
    out("Scaffolds flagged:\t", totalFlagged, "\n");
    out("Scaffolds clean:\t", totalPaths - totalFlagged, "\n");
    out("Discordant arms:\t", totalDiscordantArms, "\n");
    out("Fragmented arms:\t", totalFragmented, "\n");
}
