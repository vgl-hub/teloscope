#include <iostream>
#include <fstream>
#include <sstream> // check
#include <stdint.h> // what's this for?
#include <vector>
#include <algorithm>
#include <array> // check
#include <cmath>
#include <type_traits> // handlesBEDFile
#include <chrono> // check
#include <string_view>

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

void Trie::insertPattern(const std::string& pattern) {
    auto current = root;
    for (char ch : pattern) {
        if (current->children.find(ch) == current->children.end()) {
            current->children[ch] = std::make_shared<TrieNode>();
        }
        current = current->children[ch];
    }
    current->isEndOfWord = true;

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


std::vector<TelomereBlock> Teloscope::getBlocksRecycle(
    const std::vector<MatchInfo>& matches, 
    uint16_t mergeDist,
    std::vector<MatchInfo>& interstitialMatches,
    bool recycleToStart) {

    std::vector<TelomereBlock> telomereBlocks;
    std::vector<MatchInfo> recycledMatches;
    recycledMatches.reserve(matches.size());
    uint16_t minBlockCounts = userInput.minBlockCounts;

    // Track block indices
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

            // p/q assignment based on forward percentage
            float forwardRatio = (forwardCount * 100.0f) / blockCounts;
            if (forwardRatio > 66.0f) {
                block.blockLabel = 'p';
            } else if (forwardRatio < 33.0f) {
                block.blockLabel = 'q';
            } else {
                block.blockLabel = 'u';
            }

            telomereBlocks.push_back(block);
        } else {
            // Recycle all matches in this block
            for (size_t i = blockStartIdx; i < blockStartIdx + blockCounts; i++) {
                recycledMatches.push_back(matches[i]);
            }
        }
    };

    for (size_t i = 1; i < matches.size(); ++i) {
        uint32_t distance = matches[i].position - prevPosition;

        if (distance <= mergeDist) {
            // Continue current block
            blockCounts++;
            forwardCount += matches[i].isForward;
            canonicalCount += matches[i].isCanonical;
            
            // Update covered nucleotides - more efficient using the boolean values
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

    // Efficient move-based insertion
    if (!recycledMatches.empty()) {
        if (recycleToStart)
            interstitialMatches.insert(interstitialMatches.begin(),
                                    std::make_move_iterator(recycledMatches.begin()),
                                    std::make_move_iterator(recycledMatches.end()));
        else
            interstitialMatches.insert(interstitialMatches.end(),
                                    std::make_move_iterator(recycledMatches.begin()),
                                    std::make_move_iterator(recycledMatches.end()));
    }

    return telomereBlocks;
}


std::vector<TelomereBlock> Teloscope::getBlocks(
    std::vector<MatchInfo>& matches, 
    uint16_t mergeDist, bool needsSorting) {
    
    // Sort matches by position
    if (needsSorting) {
        std::sort(matches.begin(), matches.end(), [](const MatchInfo &a, const MatchInfo &b) {
            return a.position < b.position;
        });
    }

    std::vector<TelomereBlock> telomereBlocks;
    uint16_t minBlockCounts = userInput.minBlockCounts;

    // Initialize the first block
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

            // p/q assignment based on forward percentage
            float forwardRatio = (forwardCount * 100.0f) / blockCounts;
            if (forwardRatio > 66.6f) {
                block.blockLabel = 'p';
            } else if (forwardRatio < 33.3f) {
                block.blockLabel = 'q';
            } else {
                block.blockLabel = 'u';
            }

            telomereBlocks.push_back(block);
        }
    };

    for (size_t i = 1; i < matches.size(); ++i) {
        uint32_t distance = matches[i].position - prevPosition;
        
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

    return telomereBlocks;
}


std::vector<TelomereBlock> Teloscope::extendBlocks(std::vector<TelomereBlock> &blocks, 
    uint16_t maxBlockDist, float densityCutoff, uint32_t segmentSize, uint32_t absPos) {
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

        } else {
            // Calculate relative distances
            uint64_t relativeStart = currentBlock.start - absPos;
            uint64_t leftDist = relativeStart;
            uint64_t rightDist = segmentSize - (relativeStart + currentBlock.blockLen);
            
            // Apply overhang rules
            if ((currentBlock.blockLabel == 'p' && leftDist > rightDist) ||
                (currentBlock.blockLabel == 'q' && leftDist < rightDist)) {
                currentBlock.hasValidOr = false;
            }
            
            // Push the finalized block
            float density = static_cast<float>(currentBlock.canCovered) / currentBlock.blockLen;
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
    uint64_t rightDist = segmentSize - (relativeStart + currentBlock.blockLen);
    
    if ((currentBlock.blockLabel == 'p' && leftDist > rightDist) ||
        (currentBlock.blockLabel == 'q' && leftDist < rightDist)) {
        currentBlock.hasValidOr = false;
    }
    
    float density = static_cast<float>(currentBlock.canCovered) / currentBlock.blockLen;
    if (currentBlock.blockLen >= userInput.minBlockLen && density >= densityCutoff) {
        extendedBlocks.push_back(currentBlock);
    }

    return extendedBlocks;
}


void Teloscope::labelTerminalBlocks(
    std::vector<TelomereBlock>& blocks, uint16_t gaps,
    std::string& terminalLabel, std::string& scaffoldType) {
    // Reset isLongest for all blocks
    for (auto& block : blocks) {
        block.isLongest = false;
    }

    // Initialize return values
    terminalLabel = "";
    bool hasGaps = (gaps > 0);
    std::string gap_prefix = hasGaps ? "gapped_" : "";

    // If no blocks, set empty label and "none" type and return early
    if (blocks.empty()) {
        scaffoldType = gap_prefix + "none";
        (hasGaps ? totalGappedNone : totalNone)++; // Summary counters
        return;
    }

    // Sort blocks by coordinates for consistent ordering
    std::sort(blocks.begin(), blocks.end(),
            [](const TelomereBlock &a, const TelomereBlock &b) {
                return a.start < b.start;
            });

    // Find longest 'p' and 'q' blocks and track error status
    TelomereBlock* longest_p = nullptr;
    TelomereBlock* longest_q = nullptr;
    uint64_t max_p_coverage = 0;
    uint64_t max_q_coverage = 0;

    // Build granular label while finding longest blocks
    for (auto& block : blocks) {
        // Find longest blocks
        if (block.blockLabel == 'p' && block.canCovered > max_p_coverage) {
            longest_p = &block;
            max_p_coverage = block.canCovered;
        } 
        else if (block.blockLabel == 'q' && block.canCovered > max_q_coverage) {
            longest_q = &block;
            max_q_coverage = block.canCovered;
        }
        
        // Add to granular label (will update uppercase status later)
        char label = block.blockLabel;
        terminalLabel += label;
        
        if (!block.hasValidOr) {
            terminalLabel += '*';
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
    
    // Determine scaffold type directly from blocks info
    bool has_P = (longest_p != nullptr);
    bool has_Q = (longest_q != nullptr);
    
    // Check for discordant telomeres in longest blocks
    if ((has_P && !longest_p->hasValidOr) || (has_Q && !longest_q->hasValidOr)) {
        scaffoldType = gap_prefix + "discordant";
        (hasGaps ? totalGappedDiscordant : totalDiscordant)++; // Summary counters
        return;
    }
    
    // Complete chromosome
    if (has_P && has_Q) {
        // Check if P comes before Q in the sorted blocks
        if (longest_p->start < longest_q->start) {
            scaffoldType = gap_prefix + "t2t";
            (hasGaps ? totalGappedT2T : totalT2T)++; // Summary counters
        } else {
            // QP, Qq, Pp cases
            scaffoldType = gap_prefix + "missassembly";
            (hasGaps ? totalGappedMissassembly : totalMissassembly)++; // Summary counters
        }
        return;
    }
    
    // Check for Pp and Qq missassemblies
    if (has_P) {
        // Check for lowercase p blocks that would indicate missassembly
        bool has_p_blocks = false;
        for (const auto& block : blocks) {
            if (block.blockLabel == 'p' && &block != longest_p && block.hasValidOr) {
                has_p_blocks = true;
                break;
            }
        }
        if (has_p_blocks) {
            scaffoldType = gap_prefix + "missassembly";
            (hasGaps ? totalGappedMissassembly : totalMissassembly)++; // Summary counters
            return;
        }
    }
    
    if (has_Q) {
        // Check for lowercase q blocks that would indicate missassembly
        bool has_q_blocks = false;
        for (const auto& block : blocks) {
            if (block.blockLabel == 'q' && &block != longest_q && block.hasValidOr) {
                has_q_blocks = true;
                break;
            }
        }
        if (has_q_blocks) {
            scaffoldType = gap_prefix + "missassembly";
            (hasGaps ? totalGappedMissassembly : totalMissassembly)++; // Summary counters
            return;
        }
    }
    
    // Incomplete chromosome with a single telomere
    scaffoldType = gap_prefix + "incomplete";
    hasGaps ? totalGappedIncomplete++ : totalIncomplete++;
}


std::vector<TelomereBlock> Teloscope::filterITSBlocks(const std::vector<TelomereBlock>& interstitialBlocks) {
    std::vector<TelomereBlock> filteredBlocks;
    uint16_t minLength = 2 * userInput.patterns.front().size();

    for (const auto& block : interstitialBlocks) {
        if (block.blockLen < minLength) continue;  // Length filter
        if (block.canonicalCount == 0) continue;   // Minimal canonical check
        if ((static_cast<float>(block.canCovered) / block.blockLen) < 0.5) continue; // Majority canonical
        if (block.blockLabel == 'u' && (block.forwardCount < 2 && block.reverseCount < 2)) continue; // U check

    filteredBlocks.push_back(block);
    }

    return filteredBlocks;
}


void Teloscope::scanWindow(const std::string_view &window, uint32_t windowStart,
                           uint32_t absPos, uint32_t segmentSize,
                           WindowData& windowData, WindowData& nextOverlapData) {
    if (window.empty()) {
        return;
    }

    const unsigned short int longestPatternSize = trie.getLongestPatternSize();
    const uint32_t step = std::max<uint32_t>(1, userInput.step);
    const uint32_t overlapSize = (userInput.windowSize > step) ? (userInput.windowSize - step) : 0;
    const uint32_t terminalLimit = userInput.terminalLimit;
    const bool computeGC = userInput.outGC;
    const bool computeEntropy = userInput.outEntropy;
    int lastCanonicalPos = -1;

    // Determine starting index for Trie scanning
    uint32_t startIndex = 0;
    if (windowStart != 0 && overlapSize != 0) {
        uint32_t stepOffset = (step > longestPatternSize) ? (step - longestPatternSize) : 0;
        uint32_t overlapOffset = (overlapSize > longestPatternSize) ? (overlapSize - longestPatternSize) : 0;
        startIndex = std::min(stepOffset, overlapOffset);
    }

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

            if (i >= overlapSize || overlapSize == 0 || windowStart == 0) {
                windowData.nucleotideCounts[index]++;
            }
            if (i >= step && overlapSize != 0) {
                nextOverlapData.nucleotideCounts[index]++;
            }
        }

        // Pattern matching using Trie
        auto current = trie.getRoot();
        uint32_t scanLimit = std::min(i + longestPatternSize, static_cast<uint32_t>(window.size()));

        for (uint32_t j = i; j < scanLimit; ++j) {
            current = trie.getChild(current, window[j]);
            if (!current) break;

            if (current->isEndOfWord) {
                std::string_view pattern(window.data() + i, (j - i + 1));
                bool isTerminal = (windowStart + i <= terminalLimit ||
                                   windowStart + i >= (segmentSize > terminalLimit ? segmentSize - terminalLimit : 0));
                bool isForward = (pattern.size() >= 3 && pattern.compare(0, 3, "CCC") == 0);
                bool isCanonical = (pattern == std::string_view(userInput.canonicalFwd) ||
                                    pattern == std::string_view(userInput.canonicalRev));
                float densityGain = (static_cast<float>(pattern.size()) / window.size()) * 100.0f;
                uint32_t matchPos = absPos + windowStart + i;

                MatchInfo matchInfo;
                matchInfo.position = matchPos;
                matchInfo.isCanonical = isCanonical;
                matchInfo.isForward = isForward;
                matchInfo.matchSize = pattern.size();
                matchInfo.matchSeq = std::string(pattern);

                // Check dimers
                if (isCanonical) {
                    if (lastCanonicalPos >= 0 && (matchPos - lastCanonicalPos) <= userInput.canonicalSize) {
                        if (!windowData.hasCanDimer && (j >= overlapSize || overlapSize == 0 || windowStart == 0)) {
                            windowData.hasCanDimer = true;
                        }
                        if (!nextOverlapData.hasCanDimer && i >= step && overlapSize != 0) {
                            nextOverlapData.hasCanDimer = true;
                        }
                    }
                    lastCanonicalPos = matchPos;
                }

                if (j >= overlapSize || overlapSize == 0 || windowStart == 0) {
                    if (isCanonical) {
                        windowData.canonicalCounts++;
                        windowData.canonicalDensity += densityGain;
                    } else {
                        windowData.nonCanonicalCounts++;
                        windowData.nonCanonicalDensity += densityGain;
                    }

                    if (isForward) {
                        windowData.fwdCounts++;
                        windowData.fwdDensity += densityGain;
                    } else {
                        windowData.revCounts++;
                        windowData.revDensity += densityGain;
                    }

                    if (!isTerminal) {
                        windowData.interstitialMatches.push_back(matchInfo);
                    } else if (isForward) {
                        windowData.terminalFwdMatches.push_back(matchInfo);
                    } else {
                        windowData.terminalRevMatches.push_back(matchInfo);
                    }
                }

                if (i >= step && overlapSize != 0) {
                    if (isCanonical) {
                        nextOverlapData.canonicalCounts++;
                        nextOverlapData.canonicalDensity += densityGain;
                    } else {
                        nextOverlapData.nonCanonicalCounts++;
                        nextOverlapData.nonCanonicalDensity += densityGain;
                    }

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


SegmentData Teloscope::scanSegment(const std::string &sequence, UserInputTeloscope &userInput, uint32_t absPos) {
    SegmentData segmentData;
    const uint32_t segmentSize = sequence.size();
    if (segmentSize == 0) {
        return segmentData;
    }

    const unsigned short int longestPatternSize = trie.getLongestPatternSize();
    const uint32_t terminalLimit = userInput.terminalLimit;
    const uint32_t endBoundary = (segmentSize > terminalLimit) ? (segmentSize - terminalLimit) : 0;
    const std::string_view canonicalFwd(userInput.canonicalFwd);
    const std::string_view canonicalRev(userInput.canonicalRev);

    // Reserve common capacities to minimize reallocations
    const uint32_t estimatedMatches = std::max<uint32_t>(1, segmentSize / std::max<uint32_t>(1, userInput.canonicalSize));
    segmentData.canonicalMatches.reserve(estimatedMatches);
    segmentData.nonCanonicalMatches.reserve(estimatedMatches / 2);
    segmentData.terminalFwdMatches.reserve(estimatedMatches / 2);
    segmentData.terminalRevMatches.reserve(estimatedMatches / 2);
    segmentData.interstitialMatches.reserve(estimatedMatches);

    for (uint32_t i = 0; i < segmentSize; ++i) {
        auto node = trie.getRoot();
        uint32_t scanLimit = std::min(i + longestPatternSize, segmentSize);

        for (uint32_t j = i; j < scanLimit; ++j) {
            node = trie.getChild(node, sequence[j]);
            if (!node) break;

            if (node->isEndOfWord) {
                uint32_t len = j - i + 1;
                std::string_view pattern(sequence.data() + i, len);
                bool isForward = (len >= 3 && pattern.compare(0, 3, "CCC") == 0);
                bool isCanonical = (pattern == canonicalFwd || pattern == canonicalRev);
                uint32_t relPos = i;
                bool isTerminal = (relPos <= terminalLimit) || (relPos >= endBoundary);

                MatchInfo matchInfo;
                matchInfo.position = absPos + relPos;
                matchInfo.isCanonical = isCanonical;
                matchInfo.isForward = isForward;
                matchInfo.matchSize = len;
                matchInfo.matchSeq.assign(pattern);

                if (isCanonical) {
                    segmentData.canonicalMatches.push_back(matchInfo);
                } else if (isTerminal) {
                    segmentData.nonCanonicalMatches.push_back(matchInfo);
                }

                if (isTerminal) {
                    if (isForward) {
                        segmentData.terminalFwdMatches.push_back(matchInfo);
                    } else {
                        segmentData.terminalRevMatches.push_back(matchInfo);
                    }
                } else {
                    segmentData.interstitialMatches.push_back(matchInfo);
                }
            }
        }
    }

    // Build terminal blocks using all collected matches (distance-driven)
    const uint16_t relaxedDist = userInput.maxMatchDist;
    const uint16_t extendDist = userInput.maxBlockDist;
    const float densityCutoff = userInput.minBlockDensity;
    std::vector<TelomereBlock> fwdBlocks;
    std::vector<TelomereBlock> revBlocks;

    if (segmentData.terminalFwdMatches.size() >= 2) {
        fwdBlocks = getBlocksRecycle(segmentData.terminalFwdMatches, relaxedDist, segmentData.interstitialMatches, true);
        fwdBlocks = extendBlocks(fwdBlocks, extendDist, densityCutoff, segmentSize, absPos);
    }

    if (segmentData.terminalRevMatches.size() >= 2) {
        revBlocks = getBlocksRecycle(segmentData.terminalRevMatches, relaxedDist, segmentData.interstitialMatches, false);
        revBlocks = extendBlocks(revBlocks, extendDist, densityCutoff, segmentSize, absPos);
    }

    segmentData.terminalBlocks.reserve(fwdBlocks.size() + revBlocks.size());
    segmentData.terminalBlocks.insert(segmentData.terminalBlocks.end(), fwdBlocks.begin(), fwdBlocks.end());
    segmentData.terminalBlocks.insert(segmentData.terminalBlocks.end(), revBlocks.begin(), revBlocks.end());

    const uint16_t strictDist = 48;
    if (segmentData.interstitialMatches.size() >= 2) {
        segmentData.interstitialBlocks = getBlocks(segmentData.interstitialMatches, strictDist, true);
    }

    // Collect window metrics only when requested (non-ultra fast mode)
    if (!userInput.ultraFastMode) {
        const uint32_t step = std::max<uint32_t>(1, userInput.step);
        const uint32_t windowSize = std::max<uint32_t>(1, userInput.windowSize);
        if (segmentSize > 0) {
            uint32_t windowStart = 0;
            WindowData prevOverlapData;
            WindowData nextOverlapData;
            segmentData.windows.reserve((segmentSize / step) + 1);

            while (windowStart < segmentSize) {
                uint32_t currentWindowSize = std::min(windowSize, segmentSize - windowStart);
                std::string_view windowView(sequence.data() + windowStart, currentWindowSize);

                WindowData windowData = prevOverlapData;
                scanWindow(windowView, windowStart, absPos, segmentSize, windowData, nextOverlapData);

                if (userInput.outGC) {
                    windowData.gcContent = getGCContent(windowData.nucleotideCounts, windowView.size());
                }
                if (userInput.outEntropy) {
                    windowData.shannonEntropy = getShannonEntropy(windowData.nucleotideCounts, windowView.size());
                }

                windowData.windowStart = absPos + windowStart;
                windowData.currentWindowSize = currentWindowSize;
                segmentData.windows.emplace_back(std::move(windowData));

                prevOverlapData = nextOverlapData;
                nextOverlapData = WindowData();

                if (windowStart + step >= segmentSize) {
                    break;
                }
                windowStart += step;
            }
        }
    }

    return segmentData;
}


void Teloscope::writeBEDFile(std::ofstream& windowMetricsFile,
                            std::ofstream& canonicalMatchFile,
                            std::ofstream& noncanonicalMatchFile,
                            std::ofstream& terminalBlocksFile,
                            std::ofstream& interstitialBlocksFile) {

    // Create optimized buffers
    std::ostringstream windowMetricsBuffer;
    std::ostringstream canonicalMatchBuffer;
    std::ostringstream noncanonicalMatchBuffer;
    std::ostringstream terminalBlocksBuffer;
    std::ostringstream interstitialBlocksBuffer;

    // Write combined header (BEDgraph format)
    if (userInput.outWinRepeats || userInput.outEntropy || userInput.outGC) {
        std::string description = "";

        if (userInput.outWinRepeats) {
            description += "ForwardCounts, ReverseCounts, CanonicalCounts, NonCanonicalCounts";
        }
        if (userInput.outEntropy) {
            description += ", ShannonEntropy";
        }
        if (userInput.outGC) {
            description += ", GCContent";
        }
        windowMetricsBuffer << "track type=bedGraph name=\"Window Metrics\" description=\"" << description << "\"\n";
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
            terminalBlocksBuffer << header << "\t"
                                << block.start << "\t"
                                << blockEnd << "\t"
                                << block.blockLen << "\t"
                                << block.blockLabel << "\t"
                                << block.forwardCount << "\t"
                                << block.reverseCount << "\t"
                                << block.canonicalCount << "\t"
                                << block.nonCanonicalCount << "\t"
                                << pathSize << "\n";

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
                interstitialBlocksBuffer << header << "\t"
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
                canonicalMatchBuffer << header << "\t"
                                << match.position << "\t"
                                << (match.position + match.matchSize) << "\t"
                                << match.matchSeq << "\n";
            }

            for (const auto& match : pathData.nonCanonicalMatches) {
                noncanonicalMatchBuffer << header << "\t"
                                    << match.position << "\t"
                                    << (match.position + match.matchSize) << "\t"
                                    << match.matchSeq << "\n";
            }
        }

        // Process window metrics
        for (const auto& window : windows) {
            uint32_t windowEnd = window.windowStart + window.currentWindowSize;

            // Write window metrics if enabled
            if (userInput.outWinRepeats || userInput.outEntropy || userInput.outGC) {
                windowMetricsBuffer << header << "\t" << window.windowStart << "\t" << windowEnd;

                if (userInput.outWinRepeats) {
                    windowMetricsBuffer << "\t" << window.fwdCounts << "\t" << window.revCounts << "\t"
                                        << "\t" << window.canonicalCounts << "\t" << window.nonCanonicalCounts;
                }
                if (userInput.outEntropy) {
                    windowMetricsBuffer << "\t" << window.shannonEntropy;
                }
                if (userInput.outGC) {
                    windowMetricsBuffer << "\t" << window.gcContent;
                }
                windowMetricsBuffer << '\n';
            }
        }

        // Output path summary
        std::cout << pos + 1 << "\t" << header << "\t" 
                << longestCount << "\t"
                << (longestLabels.empty() ? "none" : longestLabels) << "\t" 
                << gaps << "\t" 
                << pathData.scaffoldType << "\t"
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

    // Single flush per file
    windowMetricsFile << windowMetricsBuffer.str();
    canonicalMatchFile << canonicalMatchBuffer.str();
    noncanonicalMatchFile << noncanonicalMatchBuffer.str();
    terminalBlocksFile << terminalBlocksBuffer.str();
    interstitialBlocksFile << interstitialBlocksBuffer.str();
}


void Teloscope::handleBEDFile() {
    lg.verbose("\nReporting window matches and metrics...");

    std::ofstream windowMetricsFile;
    std::ofstream canonicalMatchFile;
    std::ofstream noncanonicalMatchFile;
    std::ofstream terminalBlocksFile;
    std::ofstream interstitialBlocksFile;

    // Open files for writing
    if (userInput.outWinRepeats || userInput.outEntropy || userInput.outGC) {
        windowMetricsFile.open(userInput.outRoute + "/" + userInput.inSequenceName  + "_window_metrics.bedgraph");
    }

    if (userInput.outMatches) {
        canonicalMatchFile.open(userInput.outRoute + "/" + userInput.inSequenceName  + "_canonical_matches.bed");
        noncanonicalMatchFile.open(userInput.outRoute + "/" + userInput.inSequenceName  + "_noncanonical_matches.bed");
    }

    if (userInput.outITS) {
        interstitialBlocksFile.open(userInput.outRoute + "/" + userInput.inSequenceName  + "_interstitial_telomeres.bed");
    }

    terminalBlocksFile.open(userInput.outRoute + "/" + userInput.inSequenceName  + "_terminal_telomeres.bed");

    writeBEDFile(windowMetricsFile, canonicalMatchFile, noncanonicalMatchFile,
                terminalBlocksFile, interstitialBlocksFile);

    // Close all files once
    if (userInput.outWinRepeats || userInput.outEntropy || userInput.outGC) {
        windowMetricsFile.close();
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


void Teloscope::printSummary() {
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
    std::cout << "Two telomeres:\t" << totalT2T + totalGappedT2T + totalMissassembly + totalGappedMissassembly << "\n";
    std::cout << "One telomere:\t" << totalIncomplete + totalGappedIncomplete << "\n";
    std::cout << "Zero telomeres:\t" << totalNone + totalGappedNone << "\n";

    // Chromosomes by telomere completeness
    std::cout << "\n+++ Chromosome Telomere/Gap Completeness+++\n";
    std::cout << "T2T:\t" << totalT2T << "\n";
    std::cout << "Gapped T2T:\t" << totalGappedT2T << "\n";
    
    std::cout << "Missassembled:\t" << totalMissassembly << "\n";
    std::cout << "Gapped missassembled:\t" << totalGappedMissassembly << "\n";
    
    std::cout << "Incomplete:\t" << totalIncomplete << "\n";
    std::cout << "Gapped incomplete:\t" << totalGappedIncomplete << "\n";
    
    std::cout << "No telomeres:\t" << totalNone << "\n";
    std::cout << "Gapped no telomeres:\t" << totalGappedNone << "\n";
    
    std::cout << "Discordant:\t" << totalDiscordant << "\n";
    std::cout << "Gapped discordant:\t" << totalGappedDiscordant << "\n";
}