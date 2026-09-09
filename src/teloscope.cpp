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
           << " edit_distance=" << static_cast<unsigned int>(input.editDistance)
           << " ultra_fast=" << input.ultraFastMode
           << " manual_curation=" << input.manualCuration << '\n'
           << "#columns\t" << columns << '\n';
    file << header.str();
}

void writeBlockRow(std::ofstream& file, std::string_view pathName,
                   const TelomereBlock& block, uint64_t pathSize,
                   std::string_view blockType = {}) {
    const uint64_t blockEnd = block.start + block.blockLen;

    file << pathName << '\t'
         << block.start << '\t'
         << blockEnd << '\t'
         << block.blockLabel << '\t'
         << block.fwdCanCount << '\t'
         << block.revCanCount << '\t'
         << block.fwdNonCanCount << '\t'
         << block.revNonCanCount << '\t'
         << pathSize << '\t'
         << blockType << '\t'
         << block.strandLabel << '\n';
}

constexpr uint32_t minCanonicalCount = 4;

struct CoverRun {
    uint64_t start;
    uint32_t len;
};

struct Seed {
    uint64_t start;
    uint64_t end;
    uint64_t counts;
};

// per-base coverage, overlapping matches OR-ed instead of summed
void getCoverRuns(const std::vector<MatchInfo>& matches, bool canonicalOnly,
                  std::vector<CoverRun>& runs) {
    for (const MatchInfo& match : matches) {
        if (canonicalOnly && !match.isCanonical) continue;
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

// chain matches end to start, so a mixed pattern length does not shorten the gap
void getSeeds(const std::vector<MatchInfo>& matches, uint32_t matchDist,
              std::vector<Seed>& seeds) {
    for (const MatchInfo& match : matches) {
        uint64_t end = match.position + match.matchSize;
        if (!seeds.empty() && match.position <= seeds.back().end + matchDist) {
            if (end > seeds.back().end) seeds.back().end = end;
            seeds.back().counts++;
        } else {
            seeds.push_back({match.position, end, 1});
        }
    }
}

uint64_t getCoveredBases(const std::vector<CoverRun>& runs, uint64_t from, uint64_t to) {
    uint64_t covered = 0;
    for (const CoverRun& run : runs) {
        if (run.start >= to) break;
        uint64_t start = std::max(from, run.start), end = std::min(to, run.start + run.len);
        if (end > start) covered += end - start;
    }
    return covered;
}

uint64_t getGapBases(const std::vector<GapInfo>& gaps, uint64_t from, uint64_t to) {
    uint64_t total = 0;
    for (const GapInfo& gap : gaps) {
        uint64_t start = std::max(from, gap.start), end = std::min(to, gap.start + gap.length);
        if (end > start) total += end - start;
    }
    return total;
}

void tallyBlock(const std::vector<MatchInfo>& matches, TelomereBlock& block) {
    uint64_t end = block.start + block.blockLen;
    for (const MatchInfo& match : matches) {
        if (match.position >= end) break;
        if (match.position < block.start) continue;
        if (match.isCanonical) {
            if (match.isForward) block.fwdCanCount++; else block.revCanCount++;
        } else {
            if (match.isForward) block.fwdNonCanCount++; else block.revNonCanCount++;
        }
    }
}

// walk inward from the anchor and cut where the cumulative score peaks
uint64_t trimInward(const std::vector<CoverRun>& runs, const std::vector<GapInfo>& gaps,
                    uint64_t anchor, uint64_t limit, uint32_t maxGapBridge,
                    float weight, bool toRight) {
    double cumulative = 0.0, best = 0.0;
    uint64_t bestPos = anchor, prev = anchor;

    if (toRight) {
        for (const CoverRun& run : runs) {
            if (run.start < anchor) continue;
            if (run.start >= limit) break;
            uint64_t gapBases = getGapBases(gaps, prev, run.start);
            if (gapBases > maxGapBridge) break;
            cumulative -= weight * static_cast<double>((run.start - prev) - gapBases);
            cumulative += run.len;
            if (cumulative > best) { best = cumulative; bestPos = run.start + run.len; }
            prev = run.start + run.len;
        }
    } else {
        for (auto run = runs.rbegin(); run != runs.rend(); ++run) {
            uint64_t end = run->start + run->len;
            if (end > anchor) continue;
            if (end <= limit) break;
            uint64_t gapBases = getGapBases(gaps, end, prev);
            if (gapBases > maxGapBridge) break;
            cumulative -= weight * static_cast<double>((prev - end) - gapBases);
            cumulative += run->len;
            if (cumulative > best) { best = cumulative; bestPos = run->start; }
            prev = run->start;
        }
    }
    return bestPos;
}

// every locally maximal qualifying segment, so neighbouring arrays stay separate
void getMaximalSegments(const std::vector<CoverRun>& runs, const std::vector<GapInfo>& gaps,
                        uint64_t from, uint64_t to, uint32_t maxGapBridge, float weight,
                        std::vector<std::pair<uint64_t, uint64_t>>& segments) {
    double cumulative = 0.0, best = 0.0;
    uint64_t segStart = 0, segEnd = 0, prev = 0;
    bool open = false;

    auto closeSegment = [&]() {
        if (open && segEnd > segStart) segments.push_back({segStart, segEnd});
        open = false;
    };

    for (const CoverRun& run : runs) {
        uint64_t start = std::max(from, run.start);
        uint64_t end = std::min(to, run.start + run.len);
        if (run.start >= to) break;
        if (end <= start) continue;

        if (open) {
            uint64_t gapBases = getGapBases(gaps, prev, start);
            if (gapBases > maxGapBridge) {
                closeSegment();
            } else {
                cumulative -= weight * static_cast<double>((start - prev) - gapBases);
                if (cumulative < 0.0) closeSegment();
            }
        }
        if (!open) {
            open = true;
            segStart = start;
            segEnd = start;
            cumulative = 0.0;
            best = 0.0;
        }
        cumulative += static_cast<double>(end - start);
        if (cumulative > best) { best = cumulative; segEnd = end; }
        prev = end;
    }
    closeSegment();
}

} // namespace


void Teloscope::getTeloBlocks(const std::vector<MatchInfo>& matches,
                              const std::vector<GapInfo>& gapInfos, uint64_t spanSize,
                              std::vector<TelomereBlock>& terminalBlocks,
                              std::vector<TelomereBlock>& interstitialBlocks,
                              bool tipsOnly) {
    if (matches.size() < userInput.minBlockCounts) return;

    const uint32_t terminalLimit = userInput.terminalLimit;
    const uint32_t tolerance = userInput.terminalTolerance;
    const uint32_t maxGapBridge = userInput.maxBlockDist;
    const float density = userInput.minBlockDensity;
    const float weight = density / (1.0f - density);

    // first and last called base, so a leading gap cannot unanchor a real telomere
    uint64_t firstBase = 0, lastBase = spanSize;
    for (const GapInfo& gap : gapInfos) {
        if (gap.start == 0) firstBase = gap.length;
        if (gap.start + gap.length >= spanSize) lastBase = gap.start;
    }

    std::vector<CoverRun> canRuns, allRuns;
    getCoverRuns(matches, true, canRuns);
    getCoverRuns(matches, false, allRuns);

    std::vector<Seed> seeds;
    getSeeds(matches, userInput.maxMatchDist, seeds);

    auto isSeedTelomeric = [&](const Seed& seed) {
        uint64_t covered = getCoveredBases(canRuns, seed.start, seed.end);
        return seed.counts >= userInput.minBlockCounts && covered > 0 &&
               covered >= density * (seed.end - seed.start);
    };

    auto addBlock = [&](uint64_t start, uint64_t end, char arm,
                        std::vector<TelomereBlock>& out) {
        TelomereBlock block;
        block.start = start;
        block.blockLen = static_cast<uint32_t>(end - start);
        block.blockLabel = arm;
        tallyBlock(matches, block);
        const uint64_t forwardCount =
            static_cast<uint64_t>(block.fwdCanCount) + block.fwdNonCanCount;
        const uint64_t blockCounts = forwardCount +
            block.revCanCount + block.revNonCanCount;
        block.strandLabel = computeBlockLabel(forwardCount, blockCounts);
        block.hasValidOr = (block.strandLabel == arm);
        out.push_back(block);
    };

    // p arm: anchored at the first called base, trimmed inward on canonical density
    uint64_t pEnd = 0, qStart = spanSize;

    for (const Seed& seed : seeds) {
        if (seed.start > firstBase + tolerance) break;
        if (!isSeedTelomeric(seed)) continue;

        uint64_t limit = std::min<uint64_t>(spanSize, terminalLimit);
        uint64_t end = trimInward(canRuns, gapInfos, seed.start, limit,
                                  maxGapBridge, weight, true);
        uint64_t blockLen = end - seed.start;
        if (blockLen < userInput.minBlockLen) continue;
        if (getCoveredBases(canRuns, seed.start, end) < density * blockLen) continue;

        addBlock(seed.start, end, 'p', terminalBlocks);
        pEnd = end;
        break;
    }

    // q arm: the same, mirrored, and never the same physical array as the p arm
    if (spanSize > 2 * static_cast<uint64_t>(tolerance)) {
        for (auto seed = seeds.rbegin(); seed != seeds.rend(); ++seed) {
            if (seed->end + tolerance < lastBase) break;
            if (seed->start < pEnd || !isSeedTelomeric(*seed)) continue;

            uint64_t limit = (spanSize > terminalLimit) ? (spanSize - terminalLimit) : 0;
            limit = std::max(limit, pEnd);
            uint64_t start = trimInward(canRuns, gapInfos, seed->end, limit,
                                        maxGapBridge, weight, false);
            uint64_t blockLen = seed->end - start;
            if (blockLen < userInput.minBlockLen) continue;
            if (getCoveredBases(canRuns, start, seed->end) < density * blockLen) continue;

            addBlock(start, seed->end, 'q', terminalBlocks);
            qStart = start;
            break;
        }
    }

    if (tipsOnly) return;

    // interstitial blocks: variant-inclusive scoring, both edges free
    std::vector<std::pair<uint64_t, uint64_t>> segments;
    getMaximalSegments(allRuns, gapInfos, pEnd, qStart, maxGapBridge, weight, segments);

    for (const auto& segment : segments) {
        uint64_t blockLen = segment.second - segment.first;
        if (blockLen < userInput.minITSLen) continue;
        if (getCoveredBases(allRuns, segment.first, segment.second) < density * blockLen) continue;

        TelomereBlock block;
        block.start = segment.first;
        block.blockLen = static_cast<uint32_t>(blockLen);
        tallyBlock(matches, block);
        const uint64_t canonicalCount =
            static_cast<uint64_t>(block.fwdCanCount) + block.revCanCount;
        if (canonicalCount < minCanonicalCount) continue;

        const uint64_t forwardCount =
            static_cast<uint64_t>(block.fwdCanCount) + block.fwdNonCanCount;
        const uint64_t blockCounts = forwardCount + block.revCanCount + block.revNonCanCount;
        block.strandLabel = computeBlockLabel(forwardCount, blockCounts);
        // nearer end, so the column keeps one meaning across both files
        block.blockLabel = (segment.first + blockLen / 2 < spanSize / 2) ? 'p' : 'q';
        block.hasValidOr = (block.strandLabel == block.blockLabel);
        interstitialBlocks.push_back(block);
    }
}


void Teloscope::labelTerminalBlocks(
    std::vector<TelomereBlock>& blocks, uint16_t gaps,
    std::string& terminalLabel, ScaffoldType& scaffoldType,
    uint64_t pathSize, uint32_t terminalLimit) {
    (void)pathSize;
    (void)terminalLimit;

    for (auto& block : blocks) {
        block.isLongest = false;
    }

    terminalLabel = "";
    bool hasGaps = (gaps > 0);

    auto pickType = [hasGaps](ScaffoldType plain, ScaffoldType gapped) {
        return hasGaps ? gapped : plain;
    };

    if (blocks.empty()) {
        scaffoldType = pickType(ScaffoldType::NONE, ScaffoldType::GAPPED_NONE);
        return;
    }

    std::sort(blocks.begin(), blocks.end(),
            [](const TelomereBlock &a, const TelomereBlock &b) {
                return a.start < b.start;
            });

    // longest p and q, by canonical content
    TelomereBlock* longest_p = nullptr;
    TelomereBlock* longest_q = nullptr;
    uint64_t max_p_canonical_count = 0;
    uint64_t max_q_canonical_count = 0;

    for (auto& block : blocks) {
        const uint64_t canonicalCount =
            static_cast<uint64_t>(block.fwdCanCount) + block.revCanCount;
        if (block.blockLabel == 'p' && canonicalCount > max_p_canonical_count) {
            longest_p = &block;
            max_p_canonical_count = canonicalCount;
        }
        else if (block.blockLabel == 'q' && canonicalCount > max_q_canonical_count) {
            longest_q = &block;
            max_q_canonical_count = canonicalCount;
        }
    }

    if (longest_p) longest_p->isLongest = true;
    if (longest_q) longest_q->isLongest = true;

    for (auto& block : blocks) {
        terminalLabel += block.blockLabel;
        if (!block.hasValidOr) {
            terminalLabel += '*';
        }
    }

    for (size_t i = 0, j = 0; i < blocks.size(); i++) {
        if (blocks[i].isLongest) {
            terminalLabel[j] = std::toupper(terminalLabel[j]);
        }
        j++; // Move to next label position
        if (j < terminalLabel.length() && terminalLabel[j] == '*') {
            j++; // Skip asterisk if present
        }
    }

    bool has_P = (longest_p != nullptr);
    bool has_Q = (longest_q != nullptr);

    // strand balanced at an arm: real biology in some lineages, flagged not asserted
    if ((has_P && longest_p->strandLabel == 'b') || (has_Q && longest_q->strandLabel == 'b')) {
        scaffoldType = pickType(ScaffoldType::BALANCED, ScaffoldType::GAPPED_BALANCED);
        return;
    }

    // arm and strand disagree: an inverted terminal repeat
    if ((has_P && !longest_p->hasValidOr) || (has_Q && !longest_q->hasValidOr)) {
        scaffoldType = pickType(ScaffoldType::DISCORDANT, ScaffoldType::GAPPED_DISCORDANT);
        return;
    }

    // duplicated arm alongside a normal opposite arm
    for (const auto& block : blocks) {
        if (&block == longest_p || &block == longest_q) continue;
        if (block.hasValidOr) {
            scaffoldType = pickType(ScaffoldType::MISASSEMBLY, ScaffoldType::GAPPED_MISASSEMBLY);
            return;
        }
    }

    if (has_P && has_Q) {
        scaffoldType = pickType(ScaffoldType::T2T, ScaffoldType::GAPPED_T2T);
        return;
    }

    if (!has_P && !has_Q) {
        scaffoldType = pickType(ScaffoldType::NONE, ScaffoldType::GAPPED_NONE);
        return;
    }

    scaffoldType = pickType(ScaffoldType::INCOMPLETE, ScaffoldType::GAPPED_INCOMPLETE);
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

    // trie scan start index
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
        uint32_t scanLimit = std::min(i + longestPatternSize, static_cast<uint32_t>(window.size()));

        for (uint32_t j = i; j < scanLimit; ++j) {
            current = trie.getChild(current, window[j]);
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
                        if (!windowData.hasCanDimer && (alwaysMainWindow || j >= overlapSize)) {
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
                if (alwaysMainWindow || j >= overlapSize) {
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


SegmentData Teloscope::scanSegment(std::string &sequence, uint64_t absPos, bool tipsOnly) {
    SegmentData segmentData;
    uint64_t segmentSize = sequence.size();
    if (absPos + segmentSize > (1ULL << 40)) {
        std::cerr << "Error: sequence coordinate exceeds the 1.1 Tb limit.\n";
        std::exit(EXIT_FAILURE);
    }
    uint32_t terminalLimit = userInput.terminalLimit;
    unsigned short int longestPatternSize = this->trie.getLongestPatternSize();

    if (tipsOnly) {
        // ========== Fast path: terminal scan only ==========
        
        auto processRegion = [&](uint64_t start, uint64_t end) {
            for (uint64_t i = start; i < end; ++i) {
                int32_t node = trie.getRoot();
                uint64_t scanLimit = std::min(i + static_cast<uint64_t>(longestPatternSize), end);

                for (uint64_t j = i; j < scanLimit; ++j) {
                    node = trie.getChild(node, sequence[j]);
                    if (node < 0) break;

                    if (trie.isEnd(node)) {
                        uint16_t len = static_cast<uint16_t>(j - i + 1);
                        bool isForward = trie.isForward(node);
                        bool isCanonical = trie.isCanonical(node);

                        MatchInfo matchInfo{absPos + i, len, isCanonical, isForward};

                        if (isForward) segmentData.fwdCounts++;
                        else segmentData.revCounts++;
                        segmentData.allMatches.push_back(matchInfo);
                    }
                }
            }
        };
        
        if (segmentSize > 2 * terminalLimit) {
            // Process terminal regions only
            processRegion(0, terminalLimit);
            processRegion(segmentSize - terminalLimit, segmentSize);
        } else {
            // Process entire contig
            processRegion(0, segmentSize);
        }

    } else {
        // ========== Full path: window-based scan ==========
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
            ? "pos\theader\ttelomeres\tlabels\tgaps\ttype\tgranular"
            : "pos\theader\ttelomeres\tlabels\tgaps\ttype\tgranular\tits\tcanonical\twindows");

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
        std::cout << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tgranular\tits\tcanonical\twindows\n";
        reportFile << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tgranular\tits\tcanonical\twindows\n";
    } else {
        std::cout << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tgranular\n";
        reportFile << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tgranular\n";
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

        // a block reaching across an assembly gap is no longer one contig's
        auto spansGap = [&](const TelomereBlock& block) {
            for (const auto& gap : pathData.gapInfos) {
                if (gap.start < block.start + block.blockLen &&
                    gap.start + gap.length > block.start) return true;
            }
            return false;
        };

        // Terminal blocks
        std::string labels;
        for (const auto& block : pathData.terminalBlocks) {
            writeBlockRow(terminalBlocksFile, header, block, pathSize,
                          spansGap(block) ? "gapped" : "scaffold");

            // longest block only
            if (block.isLongest) {
                longestCount++;
                longestLabels += block.blockLabel;
                telomereLengths.push_back(static_cast<float>(block.blockLen));
            }
        }

        // Interstitial blocks
        if (userInput.outITS) {
            for (const auto& block : pathData.interstitialBlocks) {
                writeBlockRow(interstitialBlocksFile, header, block, pathSize,
                              spansGap(block) ? "gapped" : "contig");
            }
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
        const char* labelsStr = longestLabels.empty() ? "none" : longestLabels.c_str();

        std::cout << pos + 1 << "\t" << header << "\t"
                << longestCount << "\t" << labelsStr << "\t"
                << gaps << "\t" << typeStr << "\t"
                << pathData.terminalLabel;
        reportFile << pos + 1 << "\t" << header << "\t"
                << longestCount << "\t" << labelsStr << "\t"
                << gaps << "\t" << typeStr << "\t"
                << pathData.terminalLabel;

        totalTelomeres += longestCount;
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

    if (userInput.outITS) {
        openFile(interstitialBlocksFile, base + "_interstitial_telomeres.bed", itsBlockBuf);
    }

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

    if (userInput.outITS) {
        interstitialBlocksFile.close();
    }

    terminalBlocksFile.close();
    gapFile.close();
}


void Teloscope::computeSummaryCounts() {
    std::vector<uint64_t> scaffoldLens;
    std::vector<uint64_t> contigLens;
    scaffoldLens.reserve(allPathData.size());

    for (const auto& pathData : allPathData) {
        switch (pathData.scaffoldType) {
            case ScaffoldType::T2T:                   totalT2T++; break;
            case ScaffoldType::GAPPED_T2T:            totalGappedT2T++; break;
            case ScaffoldType::MISASSEMBLY:          totalMisassembly++; break;
            case ScaffoldType::GAPPED_MISASSEMBLY:   totalGappedMisassembly++; break;
            case ScaffoldType::INCOMPLETE:            totalIncomplete++; break;
            case ScaffoldType::GAPPED_INCOMPLETE:     totalGappedIncomplete++; break;
            case ScaffoldType::NONE:                  totalNone++; break;
            case ScaffoldType::GAPPED_NONE:           totalGappedNone++; break;
            case ScaffoldType::DISCORDANT:            totalDiscordant++; break;
            case ScaffoldType::GAPPED_DISCORDANT:     totalGappedDiscordant++; break;
            case ScaffoldType::BALANCED:              totalBalanced++; break;
            case ScaffoldType::GAPPED_BALANCED:       totalGappedBalanced++; break;
        }

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
    out("Two telomeres:\t", totalT2T + totalGappedT2T + totalMisassembly + totalGappedMisassembly, "\n");
    out("One telomere:\t", totalIncomplete + totalGappedIncomplete, "\n");
    out("Zero telomeres:\t", totalNone + totalGappedNone, "\n");

    out("\n+++ Chromosome Telomere/Gap Completeness+++\n");
    out("T2T:\t", totalT2T, "\n");
    out("Gapped T2T:\t", totalGappedT2T, "\n");

    out("Misassembled:\t", totalMisassembly, "\n");
    out("Gapped misassembled:\t", totalGappedMisassembly, "\n");

    out("Incomplete:\t", totalIncomplete, "\n");
    out("Gapped incomplete:\t", totalGappedIncomplete, "\n");

    out("No telomeres:\t", totalNone, "\n");
    out("Gapped no telomeres:\t", totalGappedNone, "\n");

    out("Discordant:\t", totalDiscordant, "\n");
    out("Gapped discordant:\t", totalGappedDiscordant, "\n");

    out("Balanced:\t", totalBalanced, "\n");
    out("Gapped balanced:\t", totalGappedBalanced, "\n");
}
