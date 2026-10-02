#include <iostream>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>

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

void writeProvenanceHeader(std::ofstream& file, const UserInputTeloscope& input,
                           std::string_view columns) {
    if (!file.is_open()) return;

    std::ostringstream header;
    header << std::boolalpha
           << "#teloscope version=" << version
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
           << " label_threshold=" << input.labelThreshold
           << " edit_distance=" << static_cast<unsigned int>(input.editDistance)
           << " ultra_fast=" << input.ultraFastMode
           << " manual_curation=" << input.manualCuration << '\n'
           << "#columns\t" << columns << '\n';
    file << header.str();
}

namespace {

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

constexpr uint32_t minCanonicalCount = 4;

enum class Orient { Fwd, Rev, All, FwdAny, RevAny };

struct Seed {
    uint64_t start;
    uint64_t end;
    uint64_t counts;
    bool isForward;
};

// per-base coverage, overlapping matches OR-ed instead of summed; a slack of -k joins one strand's matches into its chains
void getCoverRuns(const std::vector<MatchInfo>& matches, Orient orient, std::vector<CoverRun>& runs,
                  uint32_t slack = 0) {
    for (const MatchInfo& match : matches) {
        bool keep;
        switch (orient) {
            case Orient::Fwd: keep = match.isCanonical && match.isForward; break;
            case Orient::Rev: keep = match.isCanonical && !match.isForward; break;
            case Orient::FwdAny: keep = match.isForward; break;
            case Orient::RevAny: keep = !match.isForward; break;
            default:          keep = true; break;
        }
        if (!keep) continue;
        uint64_t end = match.position + match.matchSize;
        if (!runs.empty() && match.position <= runs.back().start + runs.back().len + slack) {
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

// walk the anchor's -k chain inward and cut where the cumulative score peaks, the furthest peak on a tie
uint64_t trimInward(const std::vector<CoverRun>& runs, const std::vector<CoverRun>& chains,
                    uint64_t anchor, uint64_t limit, float weight, bool toRight, uint64_t& probe) {
    double cumulative = 0.0, best = 0.0;
    uint64_t bestPos = anchor, prev = anchor;
    probe = anchor; // inner edge of the chain walked; unchanged if the anchor sits in none

    if (toRight) {
        auto chain = std::upper_bound(chains.begin(), chains.end(), anchor,
            [](uint64_t v, const CoverRun& c) { return v < c.start + c.len; });
        if (chain == chains.end() || chain->start > anchor) return bestPos;
        probe = chain->start + chain->len;
        uint64_t bound = std::min<uint64_t>(probe, limit);
        auto run = std::upper_bound(runs.begin(), runs.end(), anchor,
            [](uint64_t v, const CoverRun& r) { return v < r.start + r.len; });
        for (; run != runs.end() && run->start < bound; ++run) {
            uint64_t from = std::max(run->start, anchor);
            uint64_t end = std::min<uint64_t>(run->start + run->len, bound);
            cumulative -= weight * static_cast<double>(from - prev);
            cumulative += static_cast<double>(end - from);
            if (cumulative >= best) { best = cumulative; bestPos = end; }
            prev = end;
        }
    } else {
        auto chain = std::lower_bound(chains.begin(), chains.end(), anchor,
            [](const CoverRun& c, uint64_t v) { return c.start + c.len < v; });
        if (chain == chains.end() || chain->start >= anchor) return bestPos;
        probe = chain->start;
        uint64_t bound = std::max<uint64_t>(probe, limit);
        auto first = std::lower_bound(runs.begin(), runs.end(), anchor,
            [](const CoverRun& r, uint64_t v) { return r.start < v; });
        for (auto run = std::make_reverse_iterator(first); run != runs.rend(); ++run) {
            uint64_t end = run->start + run->len;
            if (end <= bound) break;
            uint64_t to = std::min(end, anchor);
            uint64_t start = std::max(run->start, bound);
            cumulative -= weight * static_cast<double>(prev - to);
            cumulative += static_cast<double>(to - start);
            if (cumulative >= best) { best = cumulative; bestPos = start; }
            prev = start;
        }
    }
    return bestPos;
}

// cut one seed where its score peaks, then start again past the peak, so no dense part pays for a sparse one
void getMaximalSegments(const std::vector<CoverRun>& runs, uint64_t from, uint64_t to,
                        float weight, size_t& cursor,
                        std::vector<std::pair<uint64_t, uint64_t>>& segments) {
    while (cursor < runs.size() && runs[cursor].start + runs[cursor].len <= from) ++cursor;
    size_t i = cursor;
    while (i < runs.size() && runs[i].start < to) {
        double cumulative = 0.0, best = 0.0;
        uint64_t segStart = std::max(from, runs[i].start);
        uint64_t segEnd = segStart, prev = segStart;
        size_t peak = i;
        for (size_t j = i; j < runs.size() && runs[j].start < to; ++j) {
            uint64_t start = std::max(from, runs[j].start);
            uint64_t end = std::min(to, runs[j].start + runs[j].len);
            cumulative -= weight * static_cast<double>(start - prev);
            if (cumulative < 0.0) break;
            cumulative += static_cast<double>(end - start);
            if (cumulative >= best) { best = cumulative; segEnd = end; peak = j; }
            prev = end;
        }
        if (segEnd > segStart) segments.push_back({segStart, segEnd});
        i = peak + 1;
    }
}

// junction class of each interstitial row: the nearest other-strand row within -d decides, else a same-strand row there
void assignJunctions(std::vector<TelomereBlock*>& rows, uint32_t maxBlockDist) {
    for (size_t i = 0; i < rows.size(); ++i) {
        if (rows[i]->pieces > 0) continue;
        rows[i]->junction = 's';
        const char label = rows[i]->strandLabel;
        if (label == 'b') continue;
        const uint64_t curStart = rows[i]->start, curEnd = curStart + rows[i]->blockLen;
        uint64_t bestDist = UINT64_MAX;
        bool sameStrand = false;

        // rows are start-ascending and disjoint; a mixed row ends the search on its side
        for (size_t j = i; j-- > 0; ) {
            uint64_t end = rows[j]->start + rows[j]->blockLen;
            uint64_t dist = curStart > end ? curStart - end : 0;
            char other = rows[j]->strandLabel;
            if (dist > maxBlockDist || other == 'b') break;
            if (other == label) { sameStrand = true; continue; }
            bestDist = dist;
            rows[i]->junction = (other == 'q') ? 'f' : 't';
            break;
        }
        for (size_t j = i + 1; j < rows.size(); ++j) {
            uint64_t dist = rows[j]->start > curEnd ? rows[j]->start - curEnd : 0;
            char other = rows[j]->strandLabel;
            if (dist > maxBlockDist || other == 'b') break;
            if (other == label) { sameStrand = true; continue; }
            if (dist < bestDist) {
                bestDist = dist;
                rows[i]->junction = (label == 'q') ? 'f' : 't';
            }
            break;
        }
        if (bestDist == UINT64_MAX && sameStrand) rows[i]->junction = 'r';
    }
}
} // namespace

// chr start end teloLen teloLabel closestEnd fwdCan revCan fwdNonCan revNonCan chrSize teloType
void writeBlockRow(std::ostream& file, std::string_view pathName,
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

// concordant when strandLabel matches anchorSide; complete when the read continues at least -d past the telomere
void writeReadTelomereRow(std::ostream& bedFile, const UserInputTeloscope& userInput,
                          std::string_view readName, uint64_t readLen,
                          const TelomereBlock& block, ReadTlStats& stats) {
    writeBlockRow(bedFile, readName, block, readLen, "read");
    stats.telomeresTotal++;

    if (block.strandLabel != block.anchorSide) {
        stats.telomeresDiscordant++;
        return;
    }

    uint64_t end = block.start + block.blockLen;
    uint64_t flank = (block.anchorSide == 'p') ? (readLen > end ? readLen - end : 0) : block.start;
    if (flank >= userInput.maxBlockDist) {
        stats.telomeresComplete++;
        stats.completeLengths.push_back(block.teloLen);
        stats.completeLengthSum += block.teloLen;
    } else {
        stats.telomeresReachingEnd++;
    }
}

namespace {
// R type 7 / numpy default: linear interpolation between the two closest ranks
double percentile(const std::vector<uint64_t>& sorted, double p) {
    double h = (sorted.size() - 1) * p;
    size_t lo = static_cast<size_t>(std::floor(h));
    size_t hi = std::min(lo + 1, sorted.size() - 1);
    return sorted[lo] + (h - lo) * (static_cast<double>(sorted[hi]) - sorted[lo]);
}
} // namespace

void writeReadTlReport(std::ofstream& reportFile, const UserInputTeloscope& userInput,
                       const ReadTlStats& stats) {
    writeProvenanceHeader(reportFile, userInput, "metric\tvalue");

    auto out = [&](const auto&... args) {
        std::ostringstream ss;
        (ss << ... << args);
        std::string s = ss.str();
        std::cout << s;
        reportFile << s;
    };

    out("Reads measured:\t", stats.readsMeasured, "\n");
    out("Reads kept:\t", stats.readsKept, "\n");
    out("Read telomeres:\t", stats.telomeresTotal, "\n");
    out("Complete:\t", stats.telomeresComplete, "\n");
    out("Reaching read end:\t", stats.telomeresReachingEnd, "\n");
    out("Discordant:\t", stats.telomeresDiscordant, "\n");

    if (!stats.completeLengths.empty()) {
        std::vector<uint64_t> sorted = stats.completeLengths;
        std::sort(sorted.begin(), sorted.end());
        double mean = static_cast<double>(stats.completeLengthSum) / sorted.size();

        char buf[64];
        auto fixed2 = [&](double value) {
            snprintf(buf, sizeof(buf), "%.2f", value);
            return std::string(buf);
        };

        out("Mean length:\t", fixed2(mean), "\n");
        out("Median length:\t", fixed2(percentile(sorted, 0.5)), "\n");
        out("25th percentile length:\t", fixed2(percentile(sorted, 0.25)), "\n");
        out("75th percentile length:\t", fixed2(percentile(sorted, 0.75)), "\n");
        out("90th percentile length:\t", fixed2(percentile(sorted, 0.90)), "\n");
        out("Min length:\t", sorted.front(), "\n");
        out("Max length:\t", sorted.back(), "\n");
    } else {
        out("No complete concordant telomeres for statistics.\n");
    }

    out("\nAlignment-free estimate: reads pool chromosome ends by coverage; "
        "reads broken inside a telomere look complete and pull the estimate down.\n");
}


// per contig: chain qualifying same-strand pieces inward from the anchored end, stopping at a real opposite array
TelomereBlock Teloscope::getTerminalBlocks(const std::vector<MatchInfo>& matches,
        const std::vector<CoverRun>& runsFwd, const std::vector<CoverRun>& runsRev,
        const std::vector<CoverRun>& chainsFwd, const std::vector<CoverRun>& chainsRev,
        uint64_t contigStart, uint64_t contigEnd, bool fromStart, uint64_t& outProbe) {
    TelomereBlock chain;

    const uint32_t maxBlockDist = userInput.maxBlockDist;
    const uint32_t minLen = userInput.minBlockLen;
    const uint32_t minCounts = userInput.minBlockCounts;
    const float density = userInput.minBlockDensity;
    const float weight = (density >= 1.0f) ? 1e9f : density / (1.0f - density);
    const uint32_t zone = std::min(userInput.terminalTolerance, userInput.terminalLimit);
    // a later match can still join a chain within -k or start a piece within -d
    const uint32_t margin = std::max(maxBlockDist, userInput.maxMatchDist);

    // furthest coordinate any piece or reach looked toward the interior, plus that margin
    uint64_t probeBound = fromStart ? contigStart : contigEnd;
    auto trackProbe = [&](uint64_t probe) {
        probeBound = fromStart ? std::max(probeBound, probe + margin)
                                : std::min(probeBound, probe > margin ? probe - margin : 0);
    };

    // a piece counts when it holds -c repeats of its strand and is strand-pure
    auto tallyPiece = [&](uint64_t lo, uint64_t hi, bool fwd, TelomereBlock& tally) -> bool {
        if (hi <= lo) return false;
        tally = TelomereBlock();
        tally.start = lo;
        tally.blockLen = static_cast<uint32_t>(hi - lo);
        tallyBlock(matches, tally);
        uint32_t count = fwd ? tally.fwdCanCount : tally.revCanCount;
        char label = computeStrandLabel(tally.fwdCanCount, tally.fwdCanCount + tally.revCanCount, userInput.labelThreshold);
        return count >= minCounts && label == (fwd ? 'p' : 'q');
    };

    // an array is real when unclamped pieces chained from its edge within -d sum to -l
    auto isRealArray = [&](const std::vector<CoverRun>& runs, uint64_t anchorEdge) -> bool {
        const bool isFwd = (&runs == &runsFwd);
        const std::vector<CoverRun>& chains = isFwd ? chainsFwd : chainsRev;
        const uint64_t farEdge = fromStart ? contigEnd : contigStart;
        uint64_t sum = 0, edge = anchorEdge, keptEdge = anchorEdge;
        for (bool first = true; ; first = false) {
            uint64_t probe;
            uint64_t trimmed = trimInward(runs, chains, edge, farEdge, weight, fromStart, probe);
            trackProbe(probe);
            uint64_t lo = fromStart ? edge : trimmed, hi = fromStart ? trimmed : edge;
            TelomereBlock tally;
            // a piece joins only when it outweighs the gap back to the last kept piece
            const uint64_t gap = first ? 0 : (fromStart ? lo - keptEdge : keptEdge - hi);
            if (tallyPiece(lo, hi, isFwd, tally) && static_cast<double>(hi - lo) >= weight * static_cast<double>(gap)) {
                sum += hi - lo;
                if (sum >= minLen) return true;
                keptEdge = trimmed;
                trackProbe(keptEdge);
            } else if (first) {
                return false;
            }
            // the next run of this strand past the chain just walked, while still within -d of the last kept piece
            if (fromStart) {
                auto next = std::lower_bound(runs.begin(), runs.end(), std::max(probe, edge + 1),
                    [](const CoverRun& r, uint64_t v) { return r.start < v; });
                if (next == runs.end() || next->start > keptEdge + maxBlockDist) return false;
                edge = next->start;
            } else {
                auto next = std::upper_bound(runs.begin(), runs.end(), std::min(probe, edge - 1),
                    [](uint64_t v, const CoverRun& r) { return v < r.start + r.len; });
                if (next == runs.begin()) return false;
                --next;
                if (next->start + next->len + maxBlockDist < keptEdge) return false;
                edge = next->start + next->len;
            }
        }
    };

    // trim one piece from its anchor, stopping at the first real opposite-strand array
    auto trimPiece = [&](uint64_t anchorEdge, bool fwd, uint64_t farEdge, uint64_t& chainEdge) -> uint64_t {
        const std::vector<CoverRun>& runsSame = fwd ? runsFwd : runsRev;
        const std::vector<CoverRun>& runsOpp  = fwd ? runsRev : runsFwd;
        const std::vector<CoverRun>& chains   = fwd ? chainsFwd : chainsRev;
        uint64_t probe;
        uint64_t trimmed = trimInward(runsSame, chains, anchorEdge, farEdge, weight, fromStart, probe);
        trackProbe(probe);
        chainEdge = probe;

        if (fromStart) {
            for (const CoverRun& r : runsOpp) {
                if (r.start <= anchorEdge) continue;
                if (r.start >= trimmed) break;
                if (isRealArray(runsOpp, r.start)) {
                    trimmed = trimInward(runsSame, chains, anchorEdge, r.start, weight, fromStart, probe);
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
                    trimmed = trimInward(runsSame, chains, anchorEdge, rEnd, weight, fromStart, probe);
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

    const uint64_t zoneReach = fromStart ? std::min(contigEnd, contigStart + zone)
                                          : (contigEnd > zone ? std::max(contigStart, contigEnd - zone) : contigStart);
    uint64_t reach = zoneReach;
    bool curFwd = false;
    std::vector<std::pair<uint64_t, uint64_t>> pieces;

    // a chain short of -l is dropped and the search resumes past its first piece, so a stub cannot hide a telomere
    for (size_t first = 0; first < n; ) {
        chain = TelomereBlock();
        pieces.clear();
        reach = zoneReach;
        bool haveStrand = false;
        size_t resume = n;
        uint64_t chainEdge = fromStart ? contigStart : contigEnd;
        uint64_t lastEdge = chainEdge;

        for (size_t k = first; k < n; ++k) {
            const MatchInfo& m = fromStart ? canon[k] : canon[n - 1 - k];
            uint64_t edge = fromStart ? static_cast<uint64_t>(m.position)
                                       : static_cast<uint64_t>(m.position) + m.matchSize;
            if (fromStart ? (edge > reach) : (edge < reach)) break; // a match exactly at reach is still within -d

            bool mFwd = m.isForward;
            if (haveStrand && mFwd != curFwd) {
                if (isRealArray(mFwd ? runsFwd : runsRev, edge)) break;
                continue;
            }
            // what is left of a cut chain is not this telomere: it goes to the interstitial rows
            if (haveStrand && (fromStart ? edge < chainEdge : edge > chainEdge)) continue;

            uint64_t farEdge = fromStart ? contigEnd : contigStart;
            uint64_t pieceChain;
            uint64_t trimmed = trimPiece(edge, mFwd, farEdge, pieceChain);
            uint64_t pStart = fromStart ? edge : trimmed;
            uint64_t pEnd = fromStart ? trimmed : edge;
            TelomereBlock tally;
            if (!tallyPiece(pStart, pEnd, mFwd, tally)) continue;
            // a piece joins only when it outweighs the gap back to the last piece
            if (haveStrand) {
                const uint64_t gap = fromStart ? pStart - lastEdge : lastEdge - pEnd;
                if (static_cast<double>(pEnd - pStart) < weight * static_cast<double>(gap)) continue;
            }

            pieces.push_back({pStart, pEnd});
            chain.fwdCanCount += tally.fwdCanCount;
            chain.revCanCount += tally.revCanCount;
            chain.fwdNonCanCount += tally.fwdNonCanCount;
            chain.revNonCanCount += tally.revNonCanCount;
            chain.teloLen += (pEnd - pStart);
            curFwd = mFwd;
            haveStrand = true;
            chainEdge = pieceChain;
            lastEdge = fromStart ? pEnd : pStart;
            reach = fromStart ? (pEnd + maxBlockDist) : (pStart > maxBlockDist ? pStart - maxBlockDist : 0);
            probeBound = fromStart ? std::max(probeBound, reach) : std::min(probeBound, reach);

            // past the last piece: skip everything it already consumed
            while (k + 1 < n) {
                const MatchInfo& m2 = fromStart ? canon[k + 1] : canon[n - 1 - (k + 1)];
                uint64_t edge2 = fromStart ? static_cast<uint64_t>(m2.position)
                                            : static_cast<uint64_t>(m2.position) + m2.matchSize;
                bool inside = fromStart ? (edge2 < pEnd) : (edge2 > pStart);
                if (!inside) break;
                ++k;
            }
            if (pieces.size() == 1) resume = k + 1;
        }

        if (pieces.empty() || chain.teloLen >= minLen) break;
        first = resume;
    }

    // the furthest of every piece's reach and the final anchor-search bound
    outProbe = fromStart ? std::max(probeBound, reach) : std::min(probeBound, reach);
    if (pieces.empty() || chain.teloLen < minLen) return TelomereBlock();
    uint64_t lo = pieces.front().first, hi = pieces.front().second;
    for (const auto& p : pieces) {
        lo = std::min(lo, p.first);
        hi = std::max(hi, p.second);
    }
    chain.start = lo;
    chain.blockLen = static_cast<uint32_t>(hi - lo);
    chain.pieces = static_cast<uint16_t>(std::min<size_t>(pieces.size(), UINT16_MAX));
    chain.strandLabel = curFwd ? 'p' : 'q';
    chain.anchorSide = fromStart ? 'p' : 'q';
    return chain;
}
// outside the terminal chains' spans: -k seeds with the inversion cut, each cut by getMaximalSegments
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
    size_t segCursor = 0, tallyCursor = 0, coverCursor = 0; // seeds and segments are start-ascending

    for (const Seed& seed : seeds) {
        uint64_t seedStart = std::max(seed.start, regionStart);
        uint64_t seedEnd = std::min(seed.end, regionEnd);

        std::vector<std::pair<uint64_t, uint64_t>> segments;
        if (seedEnd > seedStart) {
            getMaximalSegments(allRuns, seedStart, seedEnd, weight, segCursor, segments);
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
    }
}


void Teloscope::labelTerminalBlocks(
    std::vector<TelomereBlock>& blocks,
    ScaffoldType& scaffoldType, uint8_t& anomalyFlags) {

    anomalyFlags = 0;

    std::sort(blocks.begin(), blocks.end(),
            [](const TelomereBlock &a, const TelomereBlock &b) {
                return a.start < b.start;
            });

    TelomereBlock* armP = nullptr;
    TelomereBlock* armQ = nullptr;

    for (auto& block : blocks) {
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

    // pushes are partitioned by match start; the lookback keeps a dimer straddling the overlap edge
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
            // tiled, non-overlapping ranges extended while the terminal builder still looks that far
            uint64_t h  = buildP ? std::min<uint64_t>(t, segmentSize) : 0;
            uint64_t t0 = buildQ ? (segmentSize > t ? segmentSize - t : 0) : segmentSize;
            std::vector<MatchInfo> headMatches, tailMatches;
            if (buildP) processRegion(0, h, headMatches);
            if (buildQ) processRegion(t0, segmentSize, tailMatches);

            bool grew = true;
            while (grew) {
                grew = false;
                if (buildP && h < t0) {
                    std::vector<CoverRun> hFwd, hRev, hFwdAny, hRevAny;
                    getCoverRuns(headMatches, Orient::Fwd, hFwd);
                    getCoverRuns(headMatches, Orient::Rev, hRev);
                    getCoverRuns(headMatches, Orient::FwdAny, hFwdAny, userInput.maxMatchDist);
                    getCoverRuns(headMatches, Orient::RevAny, hRevAny, userInput.maxMatchDist);
                    uint64_t probeP;
                    getTerminalBlocks(headMatches, hFwd, hRev, hFwdAny, hRevAny, absPos, absPos + segmentSize, true, probeP);
                    if (probeP + longestPatternSize >= absPos + h) {
                        uint64_t newH = std::min<uint64_t>(h + t, t0);
                        processRegion(h, newH, headMatches);
                        h = newH;
                        grew = true;
                    }
                }
                if (buildQ && t0 > h) {
                    std::vector<CoverRun> tFwd, tRev, tFwdAny, tRevAny;
                    getCoverRuns(tailMatches, Orient::Fwd, tFwd);
                    getCoverRuns(tailMatches, Orient::Rev, tRev);
                    getCoverRuns(tailMatches, Orient::FwdAny, tFwdAny, userInput.maxMatchDist);
                    getCoverRuns(tailMatches, Orient::RevAny, tRevAny, userInput.maxMatchDist);
                    uint64_t probeQ;
                    getTerminalBlocks(tailMatches, tFwd, tRev, tFwdAny, tRevAny, absPos, absPos + segmentSize, false, probeQ);
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
    std::vector<CoverRun> runsFwd, runsRev, allRuns, chainsFwd, chainsRev;
    getCoverRuns(segmentData.allMatches, Orient::Fwd, runsFwd);
    getCoverRuns(segmentData.allMatches, Orient::Rev, runsRev);
    getCoverRuns(segmentData.allMatches, Orient::All, allRuns);
    getCoverRuns(segmentData.allMatches, Orient::FwdAny, chainsFwd, userInput.maxMatchDist);
    getCoverRuns(segmentData.allMatches, Orient::RevAny, chainsRev, userInput.maxMatchDist);

    uint64_t contigStart = absPos, contigEnd = absPos + segmentSize;
    TelomereBlock chainP, chainQ;
    uint64_t probeP, probeQ;
    if (buildP) chainP = getTerminalBlocks(segmentData.allMatches, runsFwd, runsRev, chainsFwd, chainsRev, contigStart, contigEnd, true, probeP);
    if (buildQ) chainQ = getTerminalBlocks(segmentData.allMatches, runsFwd, runsRev, chainsFwd, chainsRev, contigStart, contigEnd, false, probeQ);

    bool haveP = chainP.pieces > 0, haveQ = chainQ.pieces > 0;
    if (haveP && haveQ && chainP.strandLabel == chainQ.strandLabel &&
        chainP.start + chainP.blockLen > chainQ.start) {
        // one array reached from both ends belongs to the scaffold end, else to the end its strand points to
        if (isFirst != isLast ? isLast : chainP.strandLabel == 'q') chainP = chainQ;
        haveQ = false;
    }
    if (haveP) chainP.isScaffold = (chainP.anchorSide == 'p') ? isFirst : isLast;
    if (haveQ) chainQ.isScaffold = isLast;

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
        bool met = fastHeadActive && fastTailActive && fastHeadEnd >= fastTailStart;
        bool headIsCut = fastHeadActive && !met && fastHeadEnd != contigEnd;
        bool tailIsCut = fastTailActive && !met && fastTailStart != contigStart;
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
    assignJunctions(contigRows, userInput.maxBlockDist);

    if (haveP) segmentData.terminalBlocks.push_back(chainP);
    if (haveQ) segmentData.terminalBlocks.push_back(chainQ);
    if (userInput.outFasta) {
        if (haveP) segmentData.terminalSeqs.push_back(sequence.substr(chainP.start - absPos, chainP.blockLen));
        if (haveQ) segmentData.terminalSeqs.push_back(sequence.substr(chainQ.start - absPos, chainQ.blockLen));
    }

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
                            std::ofstream& reportFile,
                            std::ofstream& telomereFastaFile) {

    // BED and BEDGraph files carry no comment lines; the report holds the run provenance
    writeProvenanceHeader(
        reportFile, userInput,
        userInput.ultraFastMode
            ? "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly"
            : "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly\tits");

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
        std::cout << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly\tits\n";
        reportFile << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly\tits\n";
    } else {
        std::cout << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly\n";
        reportFile << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly\n";
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

        for (size_t i = 0; i < pathData.terminalBlocks.size(); ++i) {
            const auto& block = pathData.terminalBlocks[i];
            if (block.isScaffold || userInput.manualCuration) {
                writeBlockRow(terminalBlocksFile, header, block, pathSize,
                             block.isScaffold ? "scaffold" : "contig");
                if (userInput.outFasta && i < pathData.terminalSeqs.size()) {
                    telomereFastaFile << '>' << header << ':' << block.start << '-' << (block.start + block.blockLen) << '\n'
                                      << pathData.terminalSeqs[i] << '\n';
                }
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
                << anomalyStr;
        reportFile << pos + 1 << "\t" << header << "\t"
                << longestCount << "\t" << labelsStr << "\t"
                << gaps << "\t" << typeStr << "\t"
                << anomalyStr;

        totalTelomeres += longestCount;
        if (longestCount == 0) totalZeroTelomeres++;
        else if (longestCount == 1) totalOneTelomere++;
        else totalTwoTelomeres++;
        totalGaps += gaps;

        // Expand path summary
        if (!userInput.ultraFastMode) {
            std::cout << "\t" << pathData.interstitialBlocks.size();
            reportFile << "\t" << pathData.interstitialBlocks.size();
            totalITS += pathData.interstitialBlocks.size();
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
    std::vector<char> reportBuf(ioBufSize), fastaBuf(ioBufSize);

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
    std::ofstream telomereFastaFile;

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
    if (userInput.outFasta) openFile(telomereFastaFile, base + "_terminal_telomeres.fa", fastaBuf);

    writeBEDFile(windowDensityFile, windowCanonicalRatioFile,
                windowStrandRatioFile,
                windowGCFile, windowEntropyFile,
                canonicalMatchFile, noncanonicalMatchFile,
                terminalBlocksFile, interstitialBlocksFile,
                gapFile, reportFile, telomereFastaFile);

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
    if (userInput.outFasta) telomereFastaFile.close();
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

    out("\n+++ Chromosome Telomere Counts +++\n");
    out("Two telomeres:\t", totalTwoTelomeres, "\n");
    out("One telomere:\t", totalOneTelomere, "\n");
    out("Zero telomeres:\t", totalZeroTelomeres, "\n");

    // these six partition the scaffolds, which is the property worth protecting
    out("\n+++ Chromosome Telomere/Gap Completeness +++\n");
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
