#include <iostream>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits> // handlesBEDFile
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

const char* scaffoldTypeToString(ScaffoldType type) {
    switch (type) {
        case ScaffoldType::T2T:                   return "t2t";
        case ScaffoldType::GAPPED_T2T:            return "gapped_t2t";
        case ScaffoldType::MISASSEMBLY:          return "misassembly";
        case ScaffoldType::GAPPED_MISASSEMBLY:   return "gapped_misassembly";
        case ScaffoldType::INCOMPLETE:            return "incomplete";
        case ScaffoldType::GAPPED_INCOMPLETE:     return "gapped_incomplete";
        case ScaffoldType::NONE:                  return "none";
        case ScaffoldType::GAPPED_NONE:           return "gapped_none";
        case ScaffoldType::DISCORDANT:            return "discordant";
        case ScaffoldType::GAPPED_DISCORDANT:     return "gapped_discordant";
        default:                                  return "unknown";
    }
}


// Insert pattern into the flat node pool. Each char maps to index 0-3.
void Trie::insertPattern(const std::string& pattern, bool isForward, bool isCanonical) {
    int32_t current = 0; // start at root
    for (char ch : pattern) {
        int8_t idx = charToIndex(ch);
        if (idx < 0) continue; // skip non-ACGT characters
        
        if (nodes[current].children[idx] < 0) {
            // Allocate new node at end of pool
            nodes[current].children[idx] = static_cast<int32_t>(nodes.size());
            nodes.emplace_back();
        }
        current = nodes[current].children[idx];
    }
    nodes[current].isEndOfWord = true;
    nodes[current].isForward = isForward;
    nodes[current].isCanonical = isCanonical;

    if (pattern.size() > longestPatternSize) {
        longestPatternSize = pattern.size();
    }
}

Stats Teloscope::getStats(std::vector<float>& values) {
    Stats stats;
    if (values.empty()) {
        return stats;
    }

    float sum = 0.0f;
    stats.min = values[0];
    stats.max = values[0];
    for (float val : values) {
        if (val < stats.min) stats.min = val;
        if (val > stats.max) stats.max = val;
        sum += val;
    }
    stats.mean = sum / values.size();

    size_t size = values.size();
    size_t mid = size / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());

    if (size % 2 == 0) {
        float median1 = values[mid];
        float median2 = *std::max_element(values.begin(), values.begin() + mid);
        stats.median = (median1 + median2) / 2;
    } else {
        stats.median = values[mid];
    }

    return stats;
}


void Teloscope::sortBySeqPos() {
    std::sort(allPathData.begin(), allPathData.end(), [](const PathData& one, const PathData& two) {
        return one.seqPos < two.seqPos;
    });
}


char Teloscope::computeBlockLabel(uint16_t forwardCount, uint16_t blockCounts) {
    float forwardRatio = (forwardCount * 100.0f) / blockCounts;
    if (forwardRatio > 66.6f) return 'p';
    if (forwardRatio < 33.3f) return 'q';
    return 'b';
}


std::vector<TelomereBlock> Teloscope::getTeloBlocks(
    std::vector<MatchInfo>& matches, 
    uint16_t mergeDist,
    bool needsSorting,
    std::vector<MatchInfo>* recycleTarget,
    bool recycleToStart) {

    // Sort matches by position if needed
    if (needsSorting) {
        std::sort(matches.begin(), matches.end(), [](const MatchInfo &a, const MatchInfo &b) {
            return a.position < b.position;
        });
    }

    std::vector<TelomereBlock> telomereBlocks;
    std::vector<MatchInfo> recycledMatches;
    const bool doRecycle = (recycleTarget != nullptr);
    if (doRecycle) {
        recycledMatches.reserve(matches.size());
    }
    uint16_t minBlockCounts = userInput.minBlockCounts;

    // Track block indices (only needed for recycling)
    size_t blockStartIdx = 0;
    
    uint64_t blockStart = matches[0].position;
    uint64_t prevPosition = blockStart;
    uint16_t blockCounts = 1;
    uint16_t forwardCount = matches[0].isForward;
    uint16_t canonicalCount = matches[0].isCanonical;
    uint64_t totalCovered = matches[0].matchSize;
    uint64_t fwdCovered = matches[0].isForward * matches[0].matchSize;
    uint64_t canCovered = matches[0].isCanonical * matches[0].matchSize;
    uint16_t lastMatchSize = matches[0].matchSize; 

    auto finalizeBlock = [&](uint64_t endPosition) {
        if (blockCounts >= minBlockCounts && canonicalCount > 0) {
            TelomereBlock block;
            block.start = blockStart;
            block.blockLen = (endPosition - blockStart) + lastMatchSize;

            block.blockCounts = blockCounts;
            block.forwardCount = forwardCount;
            block.reverseCount = blockCounts - forwardCount;
            block.canonicalCount = canonicalCount;
            block.nonCanonicalCount = blockCounts - canonicalCount;

            block.totalCovered = totalCovered;
            block.fwdCovered = fwdCovered;
            block.canCovered = canCovered;

            block.blockLabel = computeBlockLabel(forwardCount, blockCounts);

            telomereBlocks.push_back(block);
        } else if (doRecycle) {
            // Recycle all matches in this block
            for (size_t i = blockStartIdx; i < blockStartIdx + blockCounts; i++) {
                recycledMatches.push_back(matches[i]);
            }
        }
    };

    for (size_t i = 1; i < matches.size(); ++i) {
        uint64_t distance = matches[i].position - prevPosition;

        if (distance <= mergeDist) {
            // Continue current block
            blockCounts++;
            forwardCount += matches[i].isForward;
            canonicalCount += matches[i].isCanonical;
            
            // Update covered nucleotides
            totalCovered += matches[i].matchSize;
            fwdCovered += matches[i].isForward * matches[i].matchSize;
            canCovered += matches[i].isCanonical * matches[i].matchSize;

            prevPosition = matches[i].position;
            lastMatchSize = matches[i].matchSize;

        } else {
            // Finalize current block
            finalizeBlock(prevPosition);
            
            // Start new block
            blockStartIdx = i;
            blockStart = matches[i].position;
            prevPosition = blockStart;
            blockCounts = 1;
            forwardCount = matches[i].isForward;
            canonicalCount = matches[i].isCanonical;
            
            // Reset covered nucleotides
            totalCovered = matches[i].matchSize;
            fwdCovered = matches[i].isForward * matches[i].matchSize;
            canCovered = matches[i].isCanonical * matches[i].matchSize;
            lastMatchSize = matches[i].matchSize;
        }
    }
    
    // Finalize the last block
    finalizeBlock(prevPosition);

    // Efficient move-based insertion for recycled matches
    if (doRecycle && !recycledMatches.empty()) {
        if (recycleToStart) {
            size_t oldSize = recycleTarget->size();
            recycleTarget->insert(recycleTarget->end(),
                                std::make_move_iterator(recycledMatches.begin()),
                                std::make_move_iterator(recycledMatches.end()));
            std::rotate(recycleTarget->begin(),
                        recycleTarget->begin() + oldSize,
                        recycleTarget->end());
        } else {
            recycleTarget->insert(recycleTarget->end(),
                                std::make_move_iterator(recycledMatches.begin()),
                                std::make_move_iterator(recycledMatches.end()));
        }
    }

    return telomereBlocks;
}


std::vector<TelomereBlock> Teloscope::extendBlocks(std::vector<TelomereBlock> &blocks,
    uint16_t maxBlockDist, float densityCutoff, uint64_t segmentSize, uint64_t absPos) {
    std::vector<TelomereBlock> extendedBlocks;
    if (blocks.empty()) {
        return extendedBlocks;
    }

    extendedBlocks.reserve(blocks.size());

    TelomereBlock currentBlock = blocks[0];
    currentBlock.hasValidOr = true; // Default

    for (size_t i = 1; i < blocks.size(); ++i) {
        TelomereBlock &nextBlock = blocks[i];
        uint64_t gap = nextBlock.start - (currentBlock.start + currentBlock.blockLen);

        if (gap <= maxBlockDist) {
            // Extend the current block
            currentBlock.blockLen = (nextBlock.start + nextBlock.blockLen) - currentBlock.start;

            // Update counts
            currentBlock.blockCounts += nextBlock.blockCounts;
            currentBlock.forwardCount += nextBlock.forwardCount;
            currentBlock.reverseCount += nextBlock.reverseCount;
            currentBlock.canonicalCount += nextBlock.canonicalCount;
            currentBlock.nonCanonicalCount += nextBlock.nonCanonicalCount;

            // Update coverage metrics
            currentBlock.totalCovered += nextBlock.totalCovered;
            currentBlock.fwdCovered += nextBlock.fwdCovered;
            currentBlock.canCovered += nextBlock.canCovered;

            // Recompute label from merged counts
            currentBlock.blockLabel = computeBlockLabel(currentBlock.forwardCount, currentBlock.blockCounts);

        } else {
            // Calculate relative distances (saturating to avoid underflow)
            uint64_t relativeStart = currentBlock.start - absPos;
            uint64_t leftDist = relativeStart;
            uint64_t blockEnd = relativeStart + currentBlock.blockLen;
            uint64_t rightDist = (blockEnd <= segmentSize) ? (segmentSize - blockEnd) : 0;

            // Apply overhang rules
            if ((currentBlock.blockLabel == 'p' && leftDist > rightDist) ||
                (currentBlock.blockLabel == 'q' && leftDist < rightDist)) {
                currentBlock.hasValidOr = false;
            }

            // Push the finalized block
            float density = (currentBlock.blockLen > 0)
                ? static_cast<float>(currentBlock.canCovered) / currentBlock.blockLen
                : 0.0f;
            if (currentBlock.blockLen >= userInput.minBlockLen && density >= densityCutoff) {
                extendedBlocks.push_back(currentBlock);
            }

            // Start a new block
            currentBlock = nextBlock;
            currentBlock.hasValidOr = true;  // Default is true for new blocks
        }
    }

    // Add the last block - also adjust coordinates
    uint64_t relativeStart = currentBlock.start - absPos;
    uint64_t leftDist = relativeStart;
    uint64_t blockEnd = relativeStart + currentBlock.blockLen;
    uint64_t rightDist = (blockEnd <= segmentSize) ? (segmentSize - blockEnd) : 0;

    if ((currentBlock.blockLabel == 'p' && leftDist > rightDist) ||
        (currentBlock.blockLabel == 'q' && leftDist < rightDist)) {
        currentBlock.hasValidOr = false;
    }

    float density = (currentBlock.blockLen > 0)
        ? static_cast<float>(currentBlock.canCovered) / currentBlock.blockLen
        : 0.0f;
    if (currentBlock.blockLen >= userInput.minBlockLen && density >= densityCutoff) {
        extendedBlocks.push_back(currentBlock);
    }

    return extendedBlocks;
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

    // Helper to pick gapped vs non-gapped enum variant
    auto pickType = [hasGaps](ScaffoldType plain, ScaffoldType gapped) {
        return hasGaps ? gapped : plain;
    };

    // If no blocks, set empty label and "none" type and return early
    if (blocks.empty()) {
        scaffoldType = pickType(ScaffoldType::NONE, ScaffoldType::GAPPED_NONE);
        return;
    }

    // Sort blocks by coordinates for consistent ordering
    std::sort(blocks.begin(), blocks.end(),
            [](const TelomereBlock &a, const TelomereBlock &b) {
                return a.start < b.start;
            });

    // Identify scaffold-terminal blocks (within terminalLimit of path ends)
    std::vector<TelomereBlock*> scaffoldBlocks;
    for (auto& block : blocks) {
        uint64_t blockEnd = block.start + block.blockLen;
        if (block.start < terminalLimit || blockEnd > pathSize - terminalLimit) {
            scaffoldBlocks.push_back(&block);
        }
    }

    // Build granular label from ALL blocks (for detailed annotation)
    for (auto& block : blocks) {
        terminalLabel += block.blockLabel;
        if (!block.hasValidOr) {
            terminalLabel += '*';
        }
    }

    // Find longest 'p' and 'q' blocks among SCAFFOLD-TERMINAL blocks only
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

    // Update granular label with uppercase for longest blocks
    for (size_t i = 0, j = 0; i < blocks.size(); i++) {
        if (blocks[i].isLongest) {
            terminalLabel[j] = std::toupper(terminalLabel[j]);
        }
        j++; // Move to next label position
        if (j < terminalLabel.length() && terminalLabel[j] == '*') {
            j++; // Skip asterisk if present
        }
    }

    // Determine scaffold type from scaffold-terminal blocks only
    bool has_P = (longest_p != nullptr);
    bool has_Q = (longest_q != nullptr);

    // Check for discordant telomeres in longest blocks
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

    // Check for Pp and Qq misassemblies among scaffold-terminal blocks
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

    // Incomplete chromosome with a single scaffold-terminal telomere
    scaffoldType = pickType(ScaffoldType::INCOMPLETE, ScaffoldType::GAPPED_INCOMPLETE);
}


std::vector<TelomereBlock> Teloscope::filterITSBlocks(const std::vector<TelomereBlock>& interstitialBlocks) {
    std::vector<TelomereBlock> filteredBlocks;
    uint16_t minLength = 2 * userInput.patterns.front().size();
    constexpr uint16_t minCanonicalCount = 4;

    for (const auto& block : interstitialBlocks) {
        if (block.blockLen < minLength) continue;  // Length filter
        if (block.canonicalCount < minCanonicalCount) continue; // Minimal canonical count
        if (block.blockLabel == 'b' && (block.forwardCount < 2 && block.reverseCount < 2)) continue; // Balanced check
        filteredBlocks.push_back(block);
    }

    return filteredBlocks;
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

    // Precompute densityGain lookup by match length (avoids float division per match)
    float densityByLen[256] = {};
    float invWindowSize = 100.0f / window.size();
    for (uint16_t len = 0; len <= longestPatternSize; ++len) {
        densityByLen[len] = len * invWindowSize;
    }

    // Precompute terminal status for this window (most windows are fully one or the other)
    uint64_t windowEnd = windowStart + window.size();
    uint64_t terminalEnd = (segmentSize > terminalLimit) ? (segmentSize - terminalLimit) : 0;
    bool windowFullyTerminal = (windowEnd <= terminalLimit) || (windowStart >= terminalEnd);
    bool windowFullyInterstitial = (windowStart > terminalLimit) && (windowEnd < terminalEnd);

    // Precompute overlap conditions (constant per window)
    bool alwaysMainWindow = (overlapSize == 0 || windowStart == 0);
    bool hasOverlap = (overlapSize != 0);

    // Determine starting index for Trie scanning
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

        // Pattern matching using Trie (index-based, no pointer chasing)
        int32_t current = trie.getRoot();
        uint32_t scanLimit = std::min(i + longestPatternSize, static_cast<uint32_t>(window.size()));

        for (uint32_t j = i; j < scanLimit; ++j) { // Scan positions until longest pattern
            current = trie.getChild(current, window[j]); // window[j] is a character
            if (current < 0) break; // no child → stop

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

                // Update windowData with new classification approach
                if (alwaysMainWindow || j >= overlapSize) {
                    if (isCanonical) { // C__
                        windowData.canonicalCounts++;
                        windowData.canonicalDensity += densityGain;
                        segmentData.canonicalMatches.push_back(matchInfo);
                    } else { // c__
                        windowData.nonCanonicalCounts++;
                        windowData.nonCanonicalDensity += densityGain;
                        if (isTerminal) {
                            segmentData.nonCanonicalMatches.push_back(matchInfo);
                        }
                    }

                    // Track forward/reverse metrics
                    if (isForward) {
                        windowData.fwdCounts++;
                        windowData.fwdDensity += densityGain;
                    } else {
                        windowData.revCounts++;
                        windowData.revDensity += densityGain;
                    }

                    if (!isTerminal) { // __t
                        segmentData.interstitialMatches.push_back(matchInfo);
                    } else { // __T
                        if (isForward) { // _FT
                            segmentData.terminalFwdMatches.push_back(matchInfo);
                        } else {  // _fT
                            segmentData.terminalRevMatches.push_back(matchInfo);
                        }
                    }
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

                    // Track forward/reverse metrics for overlap window
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
        // ========== FAST PATH: Direct terminal region scan ==========
        // No windows, no GC/entropy, no overlap handling - maximum speed
        
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
                            segmentData.terminalFwdMatches.push_back(matchInfo);
                        } else {
                            segmentData.terminalRevMatches.push_back(matchInfo);
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
        // ========== FULL PATH: Window-based scan with metrics ==========
        uint32_t windowSize = userInput.windowSize;
        uint32_t step = userInput.step;

        // Reserve with cap to avoid OOM on large contigs (e.g. 14GB plant chromosomes)
        constexpr uint64_t maxMatchReserve = 1000000;
        segmentData.canonicalMatches.reserve(std::min(segmentSize / 6, maxMatchReserve));
        segmentData.nonCanonicalMatches.reserve(std::min(segmentSize / 6, maxMatchReserve));
        segmentData.windows.reserve((segmentSize - windowSize) / step + 2);

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

            prevOverlapData = nextOverlapData; // Pass and reset overlap data
            nextOverlapData = WindowData(); // Reset for next iteration

            // Advance to the next window and check remaining sequence
            windowStart += step;
            if (windowStart >= segmentSize) break;

            // Prepare next window
            currentWindowSize = std::min(static_cast<uint64_t>(windowSize), segmentSize - windowStart);
            windowView = std::string_view(sequence.data() + windowStart, currentWindowSize);
        }

        segmentData.windows = std::move(windows);
    }

    // ========== Block creation from terminal matches ==========
    uint16_t mergeDist = userInput.maxMatchDist;
    uint16_t extendDist = userInput.maxBlockDist;
    float densityCutoff = userInput.minBlockDensity;
    std::vector<TelomereBlock> fwdBlocks, revBlocks;
    
    // Recycling only in full mode (nullptr disables recycling)
    std::vector<MatchInfo>* recycleTarget = tipsOnly ? nullptr : &segmentData.interstitialMatches;

    // Process forward matches
    if (segmentData.terminalFwdMatches.size() >= 2) {
        fwdBlocks = getTeloBlocks(segmentData.terminalFwdMatches, mergeDist, false, recycleTarget, true);
        fwdBlocks = extendBlocks(fwdBlocks, extendDist, densityCutoff, segmentSize, absPos);
        segmentData.terminalBlocks.insert(segmentData.terminalBlocks.end(),
                                        std::make_move_iterator(fwdBlocks.begin()),
                                        std::make_move_iterator(fwdBlocks.end()));
    }

    // Process reverse matches
    if (segmentData.terminalRevMatches.size() >= 2) {
        revBlocks = getTeloBlocks(segmentData.terminalRevMatches, mergeDist, false, recycleTarget, false);
        revBlocks = extendBlocks(revBlocks, extendDist, densityCutoff, segmentSize, absPos);
        segmentData.terminalBlocks.insert(segmentData.terminalBlocks.end(),
                                        std::make_move_iterator(revBlocks.begin()),
                                        std::make_move_iterator(revBlocks.end()));
    }

    // Create interstitial blocks (only in full mode with sufficient matches)
    if (!tipsOnly && segmentData.interstitialMatches.size() >= 2) {
        segmentData.interstitialBlocks = getTeloBlocks(segmentData.interstitialMatches, mergeDist, true);
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
                            std::ofstream& interstitialBlocksFile) {

    // Write BEDgraph headers directly to files
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

    // Console report header
    std::cout << "\n+++ Path Summary Report +++\n";
    if (!userInput.ultraFastMode) {
        std::cout << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tgranular\tits\tcanonical\twindows\n";
    } else {
        std::cout << "pos\theader\ttelomeres\tlabels\tgaps\ttype\tgranular\n";
    }

    // Processing paths
    totalPaths = allPathData.size();
    std::vector<float> telomereLengths; //  vector of floats for getStats

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

            // Determine if block is scaffold-terminal or contig-terminal
            bool isScaffoldTerminal = (block.start < userInput.terminalLimit) ||
                                      (blockEnd > pathSize - userInput.terminalLimit);
            const char* terminality = isScaffoldTerminal ? "scaffold" : "contig";

            // Default: only scaffold-terminal; --manual-curation: all
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

            // Only count and include label if this is a longest block
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

        // Output path summary
        std::cout << pos + 1 << "\t" << header << "\t" 
                << longestCount << "\t"
                << (longestLabels.empty() ? "none" : longestLabels) << "\t" 
                << gaps << "\t" 
                << scaffoldTypeToString(pathData.scaffoldType) << "\t"
                << pathData.terminalLabel;
        
        totalTelomeres += longestCount; // Update assembly summary with longest count
        totalGaps += gaps;

        // Expand path summary
        if (!userInput.ultraFastMode) {
            std::cout << "\t" 
                    << pathData.interstitialBlocks.size() << "\t"
                    << pathData.canonicalMatches.size() << "\t"
                    << windows.size();
            
            totalNWindows += windows.size(); // Update assembly summary
            totalITS += pathData.interstitialBlocks.size();
            totalCanMatches += pathData.canonicalMatches.size();
        }
        std::cout << "\n"; // Finish path summary
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

    // Buffers MUST be declared before ofstreams: C++ destroys locals in reverse
    // declaration order, so buffers outlive the streams that reference them.
    std::vector<char> densityBuf(ioBufSize), canonRatioBuf(ioBufSize), strandRatioBuf(ioBufSize);
    std::vector<char> gcBuf(ioBufSize), entropyBuf(ioBufSize);
    std::vector<char> canonMatchBuf(ioBufSize), noncanonMatchBuf(ioBufSize);
    std::vector<char> termBlockBuf(ioBufSize), itsBlockBuf(ioBufSize);

    std::ofstream windowDensityFile;
    std::ofstream windowCanonicalRatioFile;
    std::ofstream windowStrandRatioFile;
    std::ofstream windowGCFile;
    std::ofstream windowEntropyFile;
    std::ofstream canonicalMatchFile;
    std::ofstream noncanonicalMatchFile;
    std::ofstream terminalBlocksFile;
    std::ofstream interstitialBlocksFile;

    std::string base = userInput.outRoute + "/" + userInput.inSequenceName;

    // Helper to open a file with I/O buffer and check success
    // pubsetbuf with external buffers crashes on Windows (both MSVC and MinGW)
    auto openFile = [&](std::ofstream& file, const std::string& path, std::vector<char>& buf) {
#ifndef _WIN32
        file.rdbuf()->pubsetbuf(buf.data(), ioBufSize);
#else
        (void)buf;
#endif
        file.open(path);
        if (!file.is_open()) {
            fprintf(stderr, "Error: Could not open '%s' for writing.\n", path.c_str());
            exit(EXIT_FAILURE);
        }
    };

    // Open per-metric window files
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

    writeBEDFile(windowDensityFile, windowCanonicalRatioFile,
                windowStrandRatioFile,
                windowGCFile, windowEntropyFile,
                canonicalMatchFile, noncanonicalMatchFile,
                terminalBlocksFile, interstitialBlocksFile);

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


void Teloscope::printSummary() {
    computeSummaryCounts();
    std::cout << "\n+++ Assembly Summary Report +++\n";
    std::cout << "Total paths:\t" << totalPaths << "\n";
    std::cout << "Total gaps:\t" << totalGaps << "\n";
    std::cout << "Total telomeres:\t" << totalTelomeres << "\n";
    
    if (!userInput.ultraFastMode) {
        std::cout << "Total ITS blocks:\t" << totalITS << "\n";
        std::cout << "Total canonical matches:\t" << totalCanMatches << "\n";
        std::cout << "Total windows analyzed:\t" << totalNWindows << "\n";
    }
    
    // Telomere statistics
    std::cout << "\n+++ Telomere Statistics +++\n";
    if (totalTelomeres > 0) {
        std::cout << "Mean length:\t" << teloMean << "\n";
        std::cout << "Median length:\t" << teloMedian << "\n";
        std::cout << "Min length:\t" << teloMin << "\n";
        std::cout << "Max length:\t" << teloMax << "\n";
    }
    else {
        std::cout << "No telomeres found for statistics.\n";
    }

    // Chromosomes by telomere numbers (excluding discordant telomeres)
    std::cout << "\n+++ Chromosome Telomere Counts+++\n";
    std::cout << "Two telomeres:\t" << totalT2T + totalGappedT2T + totalMisassembly + totalGappedMisassembly << "\n";
    std::cout << "One telomere:\t" << totalIncomplete + totalGappedIncomplete << "\n";
    std::cout << "Zero telomeres:\t" << totalNone + totalGappedNone << "\n";

    // Chromosomes by telomere completeness
    std::cout << "\n+++ Chromosome Telomere/Gap Completeness+++\n";
    std::cout << "T2T:\t" << totalT2T << "\n";
    std::cout << "Gapped T2T:\t" << totalGappedT2T << "\n";
    
    std::cout << "Misassembled:\t" << totalMisassembly << "\n";
    std::cout << "Gapped misassembled:\t" << totalGappedMisassembly << "\n";
    
    std::cout << "Incomplete:\t" << totalIncomplete << "\n";
    std::cout << "Gapped incomplete:\t" << totalGappedIncomplete << "\n";
    
    std::cout << "No telomeres:\t" << totalNone << "\n";
    std::cout << "Gapped no telomeres:\t" << totalGappedNone << "\n";
    
    std::cout << "Discordant:\t" << totalDiscordant << "\n";
    std::cout << "Gapped discordant:\t" << totalGappedDiscordant << "\n";
}