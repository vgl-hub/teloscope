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


uint64_t Teloscope::getTerminalBlocks(
    const std::vector<MatchInfo>& matches,
    std::vector<TelomereBlock>& outBlocks,
    uint64_t segmentSize, uint64_t absPos, bool fromStart) {

    uint64_t boundary = fromStart ? absPos : (absPos + segmentSize);
    int64_t n = static_cast<int64_t>(matches.size());

    uint32_t terminalLimit = userInput.terminalLimit;
    uint16_t matchDist = userInput.maxMatchDist;
    uint16_t blockDist = userInput.maxBlockDist;
    uint16_t minBlockCounts = userInput.minBlockCounts;
    uint16_t minBlockLen = userInput.minBlockLen;
    float minBlockDensity = userInput.minBlockDensity;

    int64_t startIdx = fromStart ? 0 : n - 1;
    int64_t endIdx   = fromStart ? n : -1;
    int64_t step     = fromStart ? 1 : -1;

    auto inZone = [&](uint64_t pos) -> bool {
        uint64_t rel = pos - absPos;
        if (segmentSize <= terminalLimit) return true;
        return fromStart ? (rel < terminalLimit) : (rel >= segmentSize - terminalLimit);
    };

    // Phase 1: directional scan for sub-blocks
    std::vector<TelomereBlock> subBlocks;

    bool inBlock = false;
    uint64_t blockStart = 0, blockEnd = 0, prevPosition = 0;
    uint16_t blockCounts = 0, forwardCount = 0, canonicalCount = 0;
    uint32_t totalCovered = 0, fwdCovered = 0, canCovered = 0;

    auto startNewBlock = [&](const MatchInfo& m) {
        blockStart = m.position;
        blockEnd = m.position + m.matchSize;
        prevPosition = m.position;
        blockCounts = 1;
        forwardCount = m.isForward;
        canonicalCount = m.isCanonical;
        totalCovered = m.matchSize;
        fwdCovered = m.isForward * m.matchSize;
        canCovered = m.isCanonical * m.matchSize;
        inBlock = true;
    };

    auto finalizeSubBlock = [&]() {
        if (blockCounts >= minBlockCounts && canonicalCount > 0) {
            TelomereBlock block;
            block.start = blockStart;
            block.blockLen = static_cast<uint32_t>(blockEnd - blockStart);
            block.blockCounts = blockCounts;
            block.forwardCount = forwardCount;
            block.reverseCount = blockCounts - forwardCount;
            block.canonicalCount = canonicalCount;
            block.nonCanonicalCount = blockCounts - canonicalCount;
            block.totalCovered = totalCovered;
            block.fwdCovered = fwdCovered;
            block.canCovered = canCovered;
            subBlocks.push_back(block);
        }
        inBlock = false;
    };

    for (int64_t idx = startIdx; idx != endIdx; idx += step) {
        const MatchInfo& m = matches[idx];

        if (!inBlock) {
            if (inZone(m.position)) {
                startNewBlock(m);
            } else {
                break;
            }
        } else {
            uint64_t gap = fromStart ? (m.position - prevPosition) : (prevPosition - m.position);
            if (gap <= matchDist) {
                if (fromStart) blockEnd = m.position + m.matchSize;
                else blockStart = m.position;
                blockCounts++;
                forwardCount += m.isForward;
                canonicalCount += m.isCanonical;
                totalCovered += m.matchSize;
                fwdCovered += m.isForward * m.matchSize;
                canCovered += m.isCanonical * m.matchSize;
                prevPosition = m.position;
            } else {
                finalizeSubBlock();
                if (inZone(m.position)) {
                    startNewBlock(m);
                } else {
                    break;
                }
            }
        }
    }

    if (inBlock) finalizeSubBlock();
    if (subBlocks.empty()) return boundary;

    // Phase 2: merge sub-blocks within blockDist
    TelomereBlock current = subBlocks[0];

    auto finalizeExtended = [&]() {
        float density = (current.blockLen > 0)
            ? static_cast<float>(current.canCovered) / current.blockLen : 0.0f;

        if (current.blockLen >= minBlockLen && density >= minBlockDensity) {
            current.blockLabel = fromStart ? 'p' : 'q';

            uint64_t relStart = current.start - absPos;
            uint64_t relEnd = relStart + current.blockLen;
            uint64_t leftDist = relStart;
            uint64_t rightDist = (relEnd <= segmentSize) ? (segmentSize - relEnd) : 0;
            current.hasValidOr = fromStart ? (leftDist <= rightDist) : (leftDist >= rightDist);

            outBlocks.push_back(current);
            boundary = fromStart ? (current.start + current.blockLen) : current.start;
        }
    };

    for (size_t i = 1; i < subBlocks.size(); ++i) {
        TelomereBlock& next = subBlocks[i];
        uint64_t gap = fromStart
            ? (next.start - (current.start + current.blockLen))
            : (current.start - (next.start + next.blockLen));

        if (gap <= blockDist) {
            if (fromStart) {
                current.blockLen = (next.start + next.blockLen) - current.start;
            } else {
                current.blockLen = (current.start + current.blockLen) - next.start;
                current.start = next.start;
            }
            current.blockCounts += next.blockCounts;
            current.forwardCount += next.forwardCount;
            current.reverseCount += next.reverseCount;
            current.canonicalCount += next.canonicalCount;
            current.nonCanonicalCount += next.nonCanonicalCount;
            current.totalCovered += next.totalCovered;
            current.fwdCovered += next.fwdCovered;
            current.canCovered += next.canCovered;
        } else {
            finalizeExtended();
            current = next;
        }
    }

    finalizeExtended();
    return boundary;
}


void Teloscope::getInterstitialBlocks(
    const std::vector<MatchInfo>& allMatches,
    std::vector<TelomereBlock>& outBlocks,
    uint64_t fwdBoundary, uint64_t revBoundary) {

    uint16_t mergeDist = userInput.maxMatchDist;
    uint16_t minLength = 2 * userInput.patterns.front().size();
    constexpr uint16_t minCanonicalCount = 4;

    auto it = std::lower_bound(allMatches.begin(), allMatches.end(), fwdBoundary,
        [](const MatchInfo& m, uint64_t val) { return m.position < val; });

    if (it == allMatches.end() || it->position >= revBoundary) return;

    bool inBlock = false;
    uint64_t blockStart = 0, blockEnd = 0, prevPosition = 0;
    uint16_t blockCounts = 0, forwardCount = 0, canonicalCount = 0;
    uint32_t totalCovered = 0, fwdCovered = 0, canCovered = 0;

    auto startNewBlock = [&](const MatchInfo& m) {
        blockStart = m.position;
        blockEnd = m.position + m.matchSize;
        prevPosition = m.position;
        blockCounts = 1;
        forwardCount = m.isForward;
        canonicalCount = m.isCanonical;
        totalCovered = m.matchSize;
        fwdCovered = m.isForward * m.matchSize;
        canCovered = m.isCanonical * m.matchSize;
        inBlock = true;
    };

    auto finalizeBlock = [&]() {
        uint32_t blockLen = static_cast<uint32_t>(blockEnd - blockStart);
        char label = computeBlockLabel(forwardCount, blockCounts);

        if (blockLen >= minLength && canonicalCount >= minCanonicalCount &&
            !(label == 'b' && forwardCount < 2 && (blockCounts - forwardCount) < 2)) {

            TelomereBlock block;
            block.start = blockStart;
            block.blockLen = blockLen;
            block.blockCounts = blockCounts;
            block.forwardCount = forwardCount;
            block.reverseCount = blockCounts - forwardCount;
            block.canonicalCount = canonicalCount;
            block.nonCanonicalCount = blockCounts - canonicalCount;
            block.totalCovered = totalCovered;
            block.fwdCovered = fwdCovered;
            block.canCovered = canCovered;
            block.blockLabel = label;
            outBlocks.push_back(block);
        }
        inBlock = false;
    };

    for (; it != allMatches.end() && it->position < revBoundary; ++it) {
        const MatchInfo& m = *it;

        if (!inBlock) {
            startNewBlock(m);
        } else if (m.position - prevPosition <= mergeDist) {
            blockEnd = m.position + m.matchSize;
            blockCounts++;
            forwardCount += m.isForward;
            canonicalCount += m.isCanonical;
            totalCovered += m.matchSize;
            fwdCovered += m.isForward * m.matchSize;
            canCovered += m.isCanonical * m.matchSize;
            prevPosition = m.position;
        } else {
            finalizeBlock();
            startNewBlock(m);
        }
    }

    if (inBlock) finalizeBlock();
}


void Teloscope::labelTerminalBlocks(
    std::vector<TelomereBlock>& blocks, uint16_t gaps,
    std::string& terminalLabel, ScaffoldType& scaffoldType,
    uint64_t pathSize, uint32_t terminalLimit) {
    // Reset isLongest for all blocks
    for (auto& block : blocks) {
        block.isLongest = false;
    }

    // Initialize return values
    terminalLabel = "";
    bool hasGaps = (gaps > 0);

    auto pickType = [hasGaps](ScaffoldType plain, ScaffoldType gapped) {
        return hasGaps ? gapped : plain;
    };

    if (blocks.empty()) {
        scaffoldType = pickType(ScaffoldType::NONE, ScaffoldType::GAPPED_NONE);
        return;
    }

    // sort by coordinates
    std::sort(blocks.begin(), blocks.end(),
            [](const TelomereBlock &a, const TelomereBlock &b) {
                return a.start < b.start;
            });

    // scaffold-terminal blocks
    std::vector<TelomereBlock*> scaffoldBlocks;
    for (auto& block : blocks) {
        uint64_t blockEnd = block.start + block.blockLen;
        if (block.start < terminalLimit || blockEnd > pathSize - terminalLimit) {
            scaffoldBlocks.push_back(&block);
        }
    }

    // granular label from all blocks
    for (auto& block : blocks) {
        terminalLabel += block.blockLabel;
        if (!block.hasValidOr) {
            terminalLabel += '*';
        }
    }

    // longest p and q among scaffold-terminal blocks
    TelomereBlock* longest_p = nullptr;
    TelomereBlock* longest_q = nullptr;
    uint64_t max_p_coverage = 0;
    uint64_t max_q_coverage = 0;

    for (auto* block : scaffoldBlocks) {
        if (block->blockLabel == 'p' && block->canCovered > max_p_coverage) {
            longest_p = block;
            max_p_coverage = block->canCovered;
        }
        else if (block->blockLabel == 'q' && block->canCovered > max_q_coverage) {
            longest_q = block;
            max_q_coverage = block->canCovered;
        }
    }

    // Mark longest blocks
    if (longest_p) longest_p->isLongest = true;
    if (longest_q) longest_q->isLongest = true;

    // uppercase for longest blocks
    for (size_t i = 0, j = 0; i < blocks.size(); i++) {
        if (blocks[i].isLongest) {
            terminalLabel[j] = std::toupper(terminalLabel[j]);
        }
        j++; // Move to next label position
        if (j < terminalLabel.length() && terminalLabel[j] == '*') {
            j++; // Skip asterisk if present
        }
    }

    // scaffold type
    bool has_P = (longest_p != nullptr);
    bool has_Q = (longest_q != nullptr);

    // discordant
    if ((has_P && !longest_p->hasValidOr) || (has_Q && !longest_q->hasValidOr)) {
        scaffoldType = pickType(ScaffoldType::DISCORDANT, ScaffoldType::GAPPED_DISCORDANT);
        return;
    }

    // Complete chromosome
    if (has_P && has_Q) {
        if (longest_p->start < longest_q->start) {
            scaffoldType = pickType(ScaffoldType::T2T, ScaffoldType::GAPPED_T2T);
        } else {
            scaffoldType = pickType(ScaffoldType::MISASSEMBLY, ScaffoldType::GAPPED_MISASSEMBLY);
        }
        return;
    }

    // Pp/Qq misassemblies
    if (has_P) {
        for (const auto* block : scaffoldBlocks) {
            if (block->blockLabel == 'p' && block != longest_p && block->hasValidOr) {
                scaffoldType = pickType(ScaffoldType::MISASSEMBLY, ScaffoldType::GAPPED_MISASSEMBLY);
                return;
            }
        }
    }

    if (has_Q) {
        for (const auto* block : scaffoldBlocks) {
            if (block->blockLabel == 'q' && block != longest_q && block->hasValidOr) {
                scaffoldType = pickType(ScaffoldType::MISASSEMBLY, ScaffoldType::GAPPED_MISASSEMBLY);
                return;
            }
        }
    }

    // No scaffold-terminal telomeres found
    if (!has_P && !has_Q) {
        scaffoldType = pickType(ScaffoldType::NONE, ScaffoldType::GAPPED_NONE);
        return;
    }

    // single telomere = incomplete
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

    // density lookup by match length
    float densityByLen[256] = {};
    float invWindowSize = 100.0f / window.size();
    for (uint16_t len = 0; len <= longestPatternSize; ++len) {
        densityByLen[len] = len * invWindowSize;
    }

    // terminal status
    uint64_t windowEnd = windowStart + window.size();
    uint64_t terminalEnd = (segmentSize > terminalLimit) ? (segmentSize - terminalLimit) : 0;
    bool windowFullyTerminal = (windowEnd <= terminalLimit) || (windowStart >= terminalEnd);
    bool windowFullyInterstitial = (windowStart > terminalLimit) && (windowEnd < terminalEnd);

    // overlap conditions
    bool alwaysMainWindow = (overlapSize == 0 || windowStart == 0);
    bool hasOverlap = (overlapSize != 0);

    // trie scan start index
    uint32_t startIndex = alwaysMainWindow
                            ? 0
                            : std::min(step - longestPatternSize, overlapSize - longestPatternSize);

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
                float densityGain = densityByLen[matchLen];
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

                MatchInfo matchInfo;
                matchInfo.position = matchPos;
                matchInfo.isCanonical = isCanonical;
                matchInfo.isForward = isForward;
                matchInfo.matchSize = matchLen;
                if (needMatchSeq) {
                    matchInfo.matchSeq = std::string(window.data() + i, matchLen);
                }

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
                        windowData.canonicalDensity += densityGain;
                        segmentData.canonicalMatches.push_back(matchInfo);
                    } else {
                        windowData.nonCanonicalCounts++;
                        windowData.nonCanonicalDensity += densityGain;
                        if (isTerminal) {
                            segmentData.nonCanonicalMatches.push_back(matchInfo);
                        }
                    }

                    // route by strand
                    if (isForward) {
                        windowData.fwdCounts++;
                        windowData.fwdDensity += densityGain;
                        segmentData.fwdMatches.push_back(matchInfo);
                    } else {
                        windowData.revCounts++;
                        windowData.revDensity += densityGain;
                        segmentData.revMatches.push_back(matchInfo);
                    }
                    segmentData.allMatches.push_back(matchInfo);
                }

                // Update nextOverlapData
                if (hasOverlap && i >= step) {
                    if (isCanonical) {
                        nextOverlapData.canonicalCounts++;
                        nextOverlapData.canonicalDensity += densityGain;
                    } else {
                        nextOverlapData.nonCanonicalCounts++;
                        nextOverlapData.nonCanonicalDensity += densityGain;
                    }

                    // fwd/rev overlap metrics
                    if (isForward) {
                        nextOverlapData.fwdCounts++;
                        nextOverlapData.fwdDensity += densityGain;
                    } else {
                        nextOverlapData.revCounts++;
                        nextOverlapData.revDensity += densityGain;
                    }
                }

            }
        }
    }
}


SegmentData Teloscope::scanSegment(std::string &sequence, uint64_t absPos, bool tipsOnly) {
    SegmentData segmentData;
    uint64_t segmentSize = sequence.size();
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

                        MatchInfo matchInfo;
                        matchInfo.position = absPos + i;
                        matchInfo.isCanonical = isCanonical;
                        matchInfo.isForward = isForward;
                        matchInfo.matchSize = len;

                        if (isForward) {
                            segmentData.fwdMatches.push_back(matchInfo);
                        } else {
                            segmentData.revMatches.push_back(matchInfo);
                        }
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
        segmentData.canonicalMatches.reserve(std::min(segmentSize / 6, maxMatchReserve));
        segmentData.nonCanonicalMatches.reserve(std::min(segmentSize / 6, maxMatchReserve));
        segmentData.allMatches.reserve(std::min(segmentSize / 3, maxMatchReserve));

        if (segmentSize > windowSize) {
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
            windows.emplace_back(windowData);

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

    // ========== Block creation ==========
    uint64_t fwdBoundary = absPos;
    uint64_t revBoundary = absPos + segmentSize;

    if (segmentData.fwdMatches.size() >= 2)
        fwdBoundary = getTerminalBlocks(segmentData.fwdMatches, segmentData.terminalBlocks,
                                         segmentSize, absPos, true);
    if (segmentData.revMatches.size() >= 2)
        revBoundary = getTerminalBlocks(segmentData.revMatches, segmentData.terminalBlocks,
                                         segmentSize, absPos, false);

    if (!tipsOnly && fwdBoundary < revBoundary && segmentData.allMatches.size() >= 2)
        getInterstitialBlocks(segmentData.allMatches, segmentData.interstitialBlocks,
                              fwdBoundary, revBoundary);

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
                            std::ofstream& reportFile) {

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
        const auto& gaps = pathData.gaps;
        const auto& pathSize = pathData.pathSize;

        // Longest telomere blocks
        int longestCount = 0;
        std::string longestLabels;

        // Terminal blocks
        std::string labels;
        for (const auto& block : pathData.terminalBlocks) {
            uint64_t blockEnd = block.start + block.blockLen;

            bool isScaffoldTerminal = (block.start < userInput.terminalLimit) ||
                                      (blockEnd > pathSize - userInput.terminalLimit);
            const char* terminality = isScaffoldTerminal ? "scaffold" : "contig";

            // scaffold-terminal only, unless --manual-curation
            if (isScaffoldTerminal || userInput.manualCuration) {
                terminalBlocksFile << header << "\t"
                                    << block.start << "\t"
                                    << blockEnd << "\t"
                                    << block.blockLen << "\t"
                                    << block.blockLabel << "\t"
                                    << block.forwardCount << "\t"
                                    << block.reverseCount << "\t"
                                    << block.canonicalCount << "\t"
                                    << block.nonCanonicalCount << "\t"
                                    << pathSize << "\t"
                                    << terminality << "\n";
            }

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
                uint64_t blockEnd = block.start + block.blockLen;
                interstitialBlocksFile << header << "\t"
                                        << block.start << "\t"
                                        << blockEnd << "\t"
                                        << block.blockLen << "\t"
                                        << block.blockLabel << "\t"
                                        << block.forwardCount << "\t"
                                        << block.reverseCount << "\t"
                                        << block.canonicalCount << "\t"
                                        << block.nonCanonicalCount << "\t"
                                        << pathSize << "\n";
            }
        }

        // All canonical and terminal non-canonical matches
        if (userInput.outMatches) {
            for (const auto& match : pathData.canonicalMatches) {
                canonicalMatchFile << header << "\t"
                                << match.position << "\t"
                                << (match.position + match.matchSize) << "\t"
                                << match.matchSeq << "\n";
            }

            for (const auto& match : pathData.nonCanonicalMatches) {
                noncanonicalMatchFile << header << "\t"
                                    << match.position << "\t"
                                    << (match.position + match.matchSize) << "\t"
                                    << match.matchSeq << "\n";
            }
        }

        // Process window metrics
        for (const auto& window : windows) {
            uint64_t windowEnd = window.windowStart + window.currentWindowSize;

            if (userInput.outWinRepeats) {
                float totalDensity = window.fwdDensity + window.revDensity;
                float canonRatio = (totalDensity > 0.0f)
                    ? window.canonicalDensity / (window.canonicalDensity + window.nonCanonicalDensity)
                    : -1.0f;
                float strandRatio = (totalDensity > 0.0f)
                    ? window.fwdDensity / (window.fwdDensity + window.revDensity)
                    : -1.0f;

                windowDensityFile << header << "\t" << window.windowStart << "\t" << windowEnd
                                    << "\t" << totalDensity / 100.0f << "\n";
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
                    << pathData.canonicalMatches.size() << "\t"
                    << windows.size();
            reportFile << "\t"
                    << pathData.interstitialBlocks.size() << "\t"
                    << pathData.canonicalMatches.size() << "\t"
                    << windows.size();

            totalNWindows += windows.size();
            totalITS += pathData.interstitialBlocks.size();
            totalCanMatches += pathData.canonicalMatches.size();
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
    std::vector<char> termBlockBuf(ioBufSize), itsBlockBuf(ioBufSize);
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
    openFile(reportFile, base + "_report.tsv", reportBuf);

    writeBEDFile(windowDensityFile, windowCanonicalRatioFile,
                windowStrandRatioFile,
                windowGCFile, windowEntropyFile,
                canonicalMatchFile, noncanonicalMatchFile,
                terminalBlocksFile, interstitialBlocksFile,
                reportFile);

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
}


void Teloscope::computeSummaryCounts() {
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
        }
    }
}


void Teloscope::printSummary(std::ofstream& reportFile) {
    computeSummaryCounts();

    auto out = [&](const std::string& s) {
        std::cout << s;
        reportFile << s;
    };

    out("\n+++ Assembly Summary Report +++\n");
    out("Total paths:\t" + std::to_string(totalPaths) + "\n");
    out("Total gaps:\t" + std::to_string(totalGaps) + "\n");
    out("Total telomeres:\t" + std::to_string(totalTelomeres) + "\n");

    if (!userInput.ultraFastMode) {
        out("Total ITS blocks:\t" + std::to_string(totalITS) + "\n");
        out("Total canonical matches:\t" + std::to_string(totalCanMatches) + "\n");
        out("Total windows analyzed:\t" + std::to_string(totalNWindows) + "\n");
    }

    out("\n+++ Telomere Statistics +++\n");
    if (totalTelomeres > 0) {
        out("Mean length:\t" + std::to_string(teloMean) + "\n");
        out("Median length:\t" + std::to_string(teloMedian) + "\n");
        out("Min length:\t" + std::to_string(teloMin) + "\n");
        out("Max length:\t" + std::to_string(teloMax) + "\n");
    }
    else {
        out("No telomeres found for statistics.\n");
    }

    out("\n+++ Chromosome Telomere Counts+++\n");
    out("Two telomeres:\t" + std::to_string(totalT2T + totalGappedT2T + totalMisassembly + totalGappedMisassembly) + "\n");
    out("One telomere:\t" + std::to_string(totalIncomplete + totalGappedIncomplete) + "\n");
    out("Zero telomeres:\t" + std::to_string(totalNone + totalGappedNone) + "\n");

    out("\n+++ Chromosome Telomere/Gap Completeness+++\n");
    out("T2T:\t" + std::to_string(totalT2T) + "\n");
    out("Gapped T2T:\t" + std::to_string(totalGappedT2T) + "\n");

    out("Misassembled:\t" + std::to_string(totalMisassembly) + "\n");
    out("Gapped misassembled:\t" + std::to_string(totalGappedMisassembly) + "\n");

    out("Incomplete:\t" + std::to_string(totalIncomplete) + "\n");
    out("Gapped incomplete:\t" + std::to_string(totalGappedIncomplete) + "\n");

    out("No telomeres:\t" + std::to_string(totalNone) + "\n");
    out("Gapped no telomeres:\t" + std::to_string(totalGappedNone) + "\n");

    out("Discordant:\t" + std::to_string(totalDiscordant) + "\n");
    out("Gapped discordant:\t" + std::to_string(totalGappedDiscordant) + "\n");
}