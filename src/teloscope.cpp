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
    uint32_t counts; // matches merged into the run
};

// a telomeric stretch of one -k chain: the whole chain when it meets -y, else one of its dense segments
struct Piece {
    uint64_t start;
    uint64_t end;
    uint64_t chainStart;
    uint64_t chainEnd;
    int64_t score; // scaled: exact-repeat bases of its strand times (1 - y), minus every other base times y
    uint32_t fwdCanCount;
    uint32_t revCanCount;
    uint32_t fwdNonCanCount;
    uint32_t revNonCanCount;
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

// -y in millionths, as typed: a base inside an exact repeat scores scoreScale - density and any other base -density, so a stretch scores (its exact fraction - y) times its length, in integers
constexpr int64_t scoreScale = 1000000;
int64_t densityMillionths(float density) { return std::llround(static_cast<double>(density) * scoreScale); }

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
            runs.back().counts++;
        } else {
            runs.push_back({match.position, static_cast<uint32_t>(match.matchSize), 1});
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

// a stretch whose covered bases outweigh the rest, with its scaled score
struct ScoredSegment {
    uint64_t start;
    uint64_t end;
    int64_t score;
};

// getScoringSegments work entry: running score before and after the segment, and the last earlier entry that starts no higher
struct OpenSegment {
    uint64_t start;
    uint64_t end;
    int64_t left;
    int64_t right;
    std::ptrdiff_t lower;
};

// maximal scoring segments of the runs inside [from, to), in one pass (Ruzzo and Tompa); returns the covered bases
uint64_t getScoringSegments(const std::vector<CoverRun>& runs, size_t& cursor, uint64_t from, uint64_t to,
                            int64_t density, std::vector<OpenSegment>& open, std::vector<ScoredSegment>& segments) {
    open.clear();
    if (to <= from) return 0;
    while (cursor < runs.size() && runs[cursor].start + runs[cursor].len <= from) ++cursor;
    int64_t cumulative = 0;
    uint64_t prevEnd = from, covered = 0;
    for (size_t i = cursor; i < runs.size() && runs[i].start < to; ++i) {
        uint64_t start = std::max(from, runs[i].start);
        uint64_t end = std::min<uint64_t>(to, runs[i].start + runs[i].len);
        if (!open.empty()) cumulative -= density * static_cast<int64_t>(start - prevEnd);
        OpenSegment seg{start, end, cumulative, cumulative + (scoreScale - density) * static_cast<int64_t>(end - start), -1};
        cumulative = seg.right;
        prevEnd = end;
        covered += end - start;
        // merge leftward while an earlier segment starts no higher and ends no higher, so a tie extends the segment
        std::ptrdiff_t j = static_cast<std::ptrdiff_t>(open.size()) - 1;
        while (true) {
            while (j >= 0 && open[j].left > seg.left) j = open[j].lower;
            seg.lower = j;
            if (j < 0 || open[j].right > seg.right) break;
            seg.start = open[j].start;
            seg.left = open[j].left;
            // the search goes on below the merged entry; everything between started higher and is dropped with it
            std::ptrdiff_t next = open[j].lower;
            open.resize(static_cast<size_t>(j));
            j = next;
        }
        open.push_back(seg);
    }
    for (const OpenSegment& seg : open) segments.push_back({seg.start, seg.end, seg.right - seg.left});
    return covered;
}

// one strand's pieces: a chain of at least minCounts repeats stays whole when it meets -y, else it is cut to its dense segments
void cutChains(const std::vector<CoverRun>& runs, const std::vector<CoverRun>& chains, int64_t density,
               uint32_t minCounts, std::vector<Piece>& pieces) {
    std::vector<OpenSegment> open;
    std::vector<ScoredSegment> segments;
    size_t cursor = 0;
    for (const CoverRun& chain : chains) {
        if (chain.counts < minCounts) continue;
        const uint64_t chainEnd = chain.start + chain.len;
        segments.clear();
        uint64_t canonical = getScoringSegments(runs, cursor, chain.start, chainEnd, density, open, segments);
        if (segments.empty()) continue;
        // the whole chain, first repeat to last, variant repeats included
        const int64_t score = scoreScale * static_cast<int64_t>(canonical) - density * static_cast<int64_t>(chain.len);
        if (score >= 0) {
            pieces.push_back({chain.start, chainEnd, chain.start, chainEnd, score, 0, 0, 0, 0});
            continue;
        }
        for (const ScoredSegment& seg : segments) {
            pieces.push_back({seg.start, seg.end, chain.start, chainEnd, seg.score, 0, 0, 0, 0});
        }
    }
}

// count the repeats inside each piece by strand and kind; pieces and matches are start-ascending, so one pass serves all
void tallyPieces(const std::vector<MatchInfo>& matches, std::vector<Piece>& pieces) {
    size_t m = 0;
    for (Piece& piece : pieces) {
        while (m < matches.size() && matches[m].position < piece.start) ++m;
        for (; m < matches.size() && matches[m].position < piece.end; ++m) {
            if (matches[m].isCanonical) {
                if (matches[m].isForward) piece.fwdCanCount++; else piece.revCanCount++;
            } else {
                if (matches[m].isForward) piece.fwdNonCanCount++; else piece.revNonCanCount++;
            }
        }
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


// one strand's pieces on a contig, counted and kept when strand-pure
void Teloscope::getPieces(const std::vector<MatchInfo>& matches, const std::vector<CoverRun>& runs,
        const std::vector<CoverRun>& chains, bool isForward, std::vector<Piece>& pieces) {
    cutChains(runs, chains, densityMillionths(userInput.minBlockDensity), userInput.minBlockCounts, pieces);
    tallyPieces(matches, pieces);
    // strand-pure: more than --label-threshold of its exact repeats on its own strand
    pieces.erase(std::remove_if(pieces.begin(), pieces.end(), [&](const Piece& piece) {
        char label = computeStrandLabel(piece.fwdCanCount, piece.fwdCanCount + piece.revCanCount, userInput.labelThreshold);
        return label != (isForward ? 'p' : 'q');
    }), pieces.end());
}

// per contig end: walk each strand's pieces inward and end the telomere where their summed score peaks, the furthest peak on a tie
TelomereBlock Teloscope::getTerminalBlocks(const std::vector<Piece>& piecesFwd, const std::vector<Piece>& piecesRev,
        uint64_t contigStart, uint64_t contigEnd, bool fromStart, uint64_t& outProbe) {
    const uint32_t maxBlockDist = userInput.maxBlockDist;
    const uint32_t minLen = userInput.minBlockLen;
    const int64_t density = densityMillionths(userInput.minBlockDensity);
    const uint32_t zone = std::min(userInput.terminalTolerance, userInput.terminalLimit);
    // a later match can still join a chain within -k or start one within -d, and a chain needs --min-block-counts matches to be a piece
    const uint32_t margin = std::max(maxBlockDist, userInput.maxMatchDist) + (userInput.minBlockCounts - 1) * userInput.maxMatchDist;
    const uint64_t outerEnd = fromStart ? contigStart : contigEnd;
    const uint64_t zoneReach = fromStart ? std::min(contigEnd, contigStart + zone)
                                          : (contigEnd > zone ? std::max(contigStart, contigEnd - zone) : contigStart);

    TelomereBlock telomere;
    outProbe = zoneReach;

    // each strand's pieces in inward order: 0 forward, 1 reverse
    std::vector<const Piece*> view[2];
    const std::vector<Piece>* source[2] = {&piecesFwd, &piecesRev};
    for (int s = 0; s < 2; ++s) {
        view[s].reserve(source[s]->size());
        if (fromStart) {
            for (const Piece& piece : *source[s]) view[s].push_back(&piece);
        } else {
            for (auto it = source[s]->rbegin(); it != source[s]->rend(); ++it) view[s].push_back(&*it);
        }
    }
    if (view[0].empty() && view[1].empty()) return telomere;

    auto outer = [&](const Piece* p) { return fromStart ? p->start : p->end; };
    auto chainInner = [&](const Piece* p) { return fromStart ? p->chainEnd : p->chainStart; };
    auto before = [&](uint64_t a, uint64_t b) { return fromStart ? a < b : a > b; };
    auto further = [&](uint64_t a, uint64_t b) { return fromStart ? std::max(a, b) : std::min(a, b); };
    auto gapOf = [&](const Piece* a, const Piece* b) { return fromStart ? b->start - a->end : a->start - b->end; };
    auto spanOf = [&](const Piece* a, const Piece* b) { return fromStart ? b->end - a->start : a->end - b->start; };
    // -d is measured between the -k chains that hold two pieces
    auto joinable = [&](const Piece* a, const Piece* b) {
        return fromStart ? b->chainStart <= a->chainEnd + maxBlockDist : b->chainEnd + maxBlockDist >= a->chainStart;
    };

    // where a walk from each piece ends: within a run of joinable pieces the score is summed inward and its last peak kept
    auto walkEnds = [&](const std::vector<const Piece*>& v, const std::vector<uint8_t>& cut,
                        std::vector<uint32_t>& best, std::vector<uint32_t>& last) {
        const size_t n = v.size();
        best.assign(n, 0);
        last.assign(n, 0);
        std::vector<int64_t> sum(n);
        for (size_t s = 0; s < n; ) {
            size_t e = s;
            sum[s] = v[s]->score;
            while (e + 1 < n && !cut[e + 1]) {
                sum[e + 1] = sum[e] - density * static_cast<int64_t>(gapOf(v[e], v[e + 1])) + v[e + 1]->score;
                ++e;
            }
            size_t top = e;
            for (size_t i = e + 1; i-- > s; ) {
                if (sum[i] > sum[top]) top = i;
                best[i] = static_cast<uint32_t>(top);
                last[i] = static_cast<uint32_t>(e);
            }
            s = e + 1;
        }
    };

    // an array is real when a walk over its own strand alone, from that piece, spans -l or more
    std::vector<uint8_t> cut[2], real[2];
    std::vector<uint32_t> best[2], last[2];
    for (int s = 0; s < 2; ++s) {
        const std::vector<const Piece*>& v = view[s];
        cut[s].assign(v.size(), 1);
        for (size_t i = 1; i < v.size(); ++i) cut[s][i] = !joinable(v[i - 1], v[i]);
        walkEnds(v, cut[s], best[s], last[s]);
        real[s].resize(v.size());
        for (size_t i = 0; i < v.size(); ++i) real[s][i] = spanOf(v[i], v[best[s][i]]) >= minLen;
    }

    // a telomere also stops at a real array of the other strand; looked is how far inward each walk depended on the sequence
    std::vector<uint32_t> walkEnd[2];
    std::vector<uint64_t> looked[2];
    for (int s = 0; s < 2; ++s) {
        const std::vector<const Piece*>& v = view[s];
        const std::vector<const Piece*>& w = view[1 - s];
        const size_t n = v.size();
        std::vector<uint8_t> stop(cut[s]);
        std::vector<uint64_t> seen(n, outerEnd);
        size_t j = 0;
        for (size_t i = 1; i < n; ++i) {
            while (j < w.size() && before(outer(w[j]), outer(v[i]))) {
                if (!cut[s][i] && !before(outer(w[j]), outer(v[i - 1]))) {
                    if (real[1 - s][j]) stop[i] = 1;
                    seen[i] = further(seen[i], chainInner(w[last[1 - s][j]]));
                }
                ++j;
            }
        }
        std::vector<uint32_t> segLast;
        walkEnds(v, stop, walkEnd[s], segLast);
        looked[s].assign(n, outerEnd);
        uint64_t behind = outerEnd;
        for (size_t i = n; i-- > 0; ) {
            if (i == segLast[i]) behind = (i + 1 < n && !cut[s][i + 1]) ? seen[i + 1] : outerEnd;
            else behind = further(behind, seen[i + 1]);
            looked[s][i] = further(chainInner(v[segLast[i]]), behind);
        }
    }

    // a chain that starts in the zone must be scanned whole before its density is judged, wherever its pieces begin
    uint64_t probe = outerEnd;
    for (int s = 0; s < 2; ++s) {
        for (const Piece* p : view[s]) {
            if (before(zoneReach, fromStart ? p->chainStart : p->chainEnd)) break;
            probe = further(probe, chainInner(p));
        }
    }

    // the first piece of either strand in the start zone whose walk spans -l is the telomere; a shorter walk is skipped, so a stub cannot hide a telomere
    size_t next[2] = {0, 0};
    while (true) {
        int s = -1;
        for (int c = 0; c < 2; ++c) {
            if (next[c] >= view[c].size()) continue;
            if (s < 0 || before(outer(view[c][next[c]]), outer(view[s][next[s]]))) s = c;
        }
        if (s < 0) break;
        const size_t i = next[s]++;
        const Piece* anchor = view[s][i];
        if (before(zoneReach, outer(anchor))) break;
        probe = further(probe, looked[s][i]);
        const size_t e = walkEnd[s][i];
        if (spanOf(anchor, view[s][e]) < minLen) continue;

        uint64_t lo = anchor->start, hi = anchor->end;
        for (size_t k = i; k <= e; ++k) {
            const Piece* p = view[s][k];
            lo = std::min(lo, p->start);
            hi = std::max(hi, p->end);
            telomere.teloLen += p->end - p->start;
            telomere.fwdCanCount += p->fwdCanCount;
            telomere.revCanCount += p->revCanCount;
            telomere.fwdNonCanCount += p->fwdNonCanCount;
            telomere.revNonCanCount += p->revNonCanCount;
        }
        telomere.start = lo;
        telomere.blockLen = static_cast<uint32_t>(hi - lo);
        telomere.pieces = static_cast<uint16_t>(std::min<size_t>(e - i + 1, UINT16_MAX));
        telomere.strandLabel = (s == 0) ? 'p' : 'q';
        telomere.anchorSide = fromStart ? 'p' : 'q';
        break;
    }

    // the scan must reach past everything a walk depended on, and past the start zone
    if (probe != outerEnd) probe = fromStart ? probe + margin : (probe > margin ? probe - margin : 0);
    outProbe = further(probe, zoneReach);
    return telomere;
}

// outside the telomeres: -k seeds with the inversion cut, each cut into its maximal scoring segments on the coverage of all matches
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

    const int64_t density = densityMillionths(userInput.minBlockDensity);
    size_t segCursor = 0, tallyCursor = 0; // seeds and segments are start-ascending
    std::vector<OpenSegment> open;
    std::vector<ScoredSegment> segments;

    for (const Seed& seed : seeds) {
        uint64_t seedStart = std::max(seed.start, regionStart);
        uint64_t seedEnd = std::min(seed.end, regionEnd);

        segments.clear();
        if (seedEnd > seedStart) {
            getScoringSegments(allRuns, segCursor, seedStart, seedEnd, density, open, segments);
        }

        for (const ScoredSegment& seg : segments) {
            TelomereBlock block;
            block.start = seg.start;
            block.blockLen = static_cast<uint32_t>(seg.end - seg.start);
            tallyBlock(matches, block, tallyCursor);
            const uint64_t canonicalCount =
                static_cast<uint64_t>(block.fwdCanCount) + block.revCanCount;
            if (canonicalCount < minCanonicalCount) continue;

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
            // tiled, non-overlapping ranges, widened while the terminal builder still looks that far
            uint64_t h  = buildP ? std::min<uint64_t>(t, segmentSize) : 0;
            uint64_t t0 = buildQ ? (segmentSize > t ? segmentSize - t : 0) : segmentSize;
            std::vector<MatchInfo> headMatches, tailMatches;
            if (buildP) processRegion(0, h, headMatches);
            if (buildQ) processRegion(t0, segmentSize, tailMatches);

            // how far the terminal builder looks from one end of the matches read so far
            auto probeEnd = [&](const std::vector<MatchInfo>& ms, bool fromStart) {
                std::vector<CoverRun> fwd, rev, fwdAny, revAny;
                getCoverRuns(ms, Orient::Fwd, fwd);
                getCoverRuns(ms, Orient::Rev, rev);
                getCoverRuns(ms, Orient::FwdAny, fwdAny, userInput.maxMatchDist);
                getCoverRuns(ms, Orient::RevAny, revAny, userInput.maxMatchDist);
                std::vector<Piece> piecesFwd, piecesRev;
                getPieces(ms, fwd, fwdAny, true, piecesFwd);
                getPieces(ms, rev, revAny, false, piecesRev);
                uint64_t probe;
                getTerminalBlocks(piecesFwd, piecesRev, absPos, absPos + segmentSize, fromStart, probe);
                return probe;
            };

            // each widening doubles the last, so all rebuilds together cost a constant times the final range
            uint64_t probeP = absPos, probeQ = absPos + segmentSize, stepP = t, stepQ = t;
            bool grew = true;
            while (grew) {
                grew = false;
                if (buildP && h < t0) {
                    probeP = probeEnd(headMatches, true);
                    if (probeP + longestPatternSize >= absPos + h) {
                        uint64_t newH = std::min<uint64_t>(h + stepP, t0);
                        processRegion(h, newH, headMatches);
                        h = newH;
                        stepP *= 2;
                        grew = true;
                    }
                }
                if (buildQ && t0 > h) {
                    probeQ = probeEnd(tailMatches, false);
                    if (probeQ <= absPos + t0 + longestPatternSize) {
                        uint64_t newT0 = (t0 > stepQ) ? std::max<uint64_t>(t0 - stepQ, h) : h;
                        std::vector<MatchInfo> extra;
                        processRegion(newT0, t0, extra);
                        extra.insert(extra.end(), tailMatches.begin(), tailMatches.end());
                        tailMatches.swap(extra);
                        t0 = newT0;
                        stepQ *= 2;
                        grew = true;
                    }
                }
            }

            // the ranges met: the last probes are stale, so take them from the whole contig
            if (h >= t0) {
                std::vector<MatchInfo> all(headMatches);
                all.insert(all.end(), tailMatches.begin(), tailMatches.end());
                if (buildP) probeP = probeEnd(all, true);
                if (buildQ) probeQ = probeEnd(all, false);
            }
            // doubling can overshoot: keep the smallest multiple of -t each end needed, so the scan matches a growth by -t
            uint64_t hNeed = h, t0Need = t0;
            if (buildP) {
                uint64_t reach = probeP + longestPatternSize - absPos;
                hNeed = std::min<uint64_t>((reach / t + 1) * t, segmentSize);
            }
            if (buildQ) {
                uint64_t limit = probeQ > absPos + longestPatternSize ? probeQ - absPos - longestPatternSize : 0;
                uint64_t back = ((segmentSize - limit) / t + 1) * t;
                t0Need = back >= segmentSize ? 0 : segmentSize - back;
            }
            if (hNeed < t0Need) {
                if (buildP && hNeed < h) {
                    auto cutAt = std::lower_bound(headMatches.begin(), headMatches.end(), absPos + hNeed,
                        [](const MatchInfo& m, uint64_t v) { return m.position < v; });
                    headMatches.erase(cutAt, headMatches.end());
                    h = hNeed;
                }
                if (buildQ && t0Need > t0) {
                    auto cutAt = std::lower_bound(tailMatches.begin(), tailMatches.end(), absPos + t0Need,
                        [](const MatchInfo& m, uint64_t v) { return m.position < v; });
                    tailMatches.erase(tailMatches.begin(), cutAt);
                    t0 = t0Need;
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
    std::vector<Piece> piecesFwd, piecesRev;
    getPieces(segmentData.allMatches, runsFwd, chainsFwd, true, piecesFwd);
    getPieces(segmentData.allMatches, runsRev, chainsRev, false, piecesRev);

    uint64_t contigStart = absPos, contigEnd = absPos + segmentSize;
    TelomereBlock chainP, chainQ;
    uint64_t probeP, probeQ;
    if (buildP) chainP = getTerminalBlocks(piecesFwd, piecesRev, contigStart, contigEnd, true, probeP);
    if (buildQ) chainQ = getTerminalBlocks(piecesFwd, piecesRev, contigStart, contigEnd, false, probeQ);

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
        // fast mode: drop interstitial rows whose -k chain a wider window could still extend, only at a real cut
        const uint64_t margin = userInput.maxBlockDist + userInput.maxMatchDist + longestPatternSize;
        bool met = fastHeadActive && fastTailActive && fastHeadEnd >= fastTailStart;
        bool headIsCut = fastHeadActive && !met && fastHeadEnd != contigEnd;
        bool tailIsCut = fastTailActive && !met && fastTailStart != contigStart;
        std::vector<CoverRun> allChains;
        getCoverRuns(segmentData.allMatches, Orient::All, allChains, userInput.maxMatchDist);
        auto& rows = segmentData.interstitialBlocks;
        rows.erase(std::remove_if(rows.begin(), rows.end(), [&](const TelomereBlock& b) {
            // a row is cut from a seed, and a seed truncated at the cut can shift every row cut from it
            auto chain = std::upper_bound(allChains.begin(), allChains.end(), b.start,
                [](uint64_t v, const CoverRun& c) { return v < c.start + c.len; });
            uint64_t from = b.start, to = b.start + b.blockLen;
            if (chain != allChains.end()) {
                from = std::min<uint64_t>(from, chain->start);
                to = std::max<uint64_t>(to, chain->start + chain->len);
            }
            bool nearHead = headIsCut && b.start < fastHeadEnd && (to + margin >= fastHeadEnd);
            bool nearTail = tailIsCut && b.start >= fastTailStart && (from <= fastTailStart + margin);
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
