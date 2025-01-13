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


std::vector<TelomereBlock> Teloscope::getTelomereBlocks(const std::vector<MatchInfo>& inputMatches, uint16_t mergeDist) {
    std::vector<TelomereBlock> telomereBlocks;
    uint16_t patternSize = userInput.canonicalSize;
    uint16_t minBlockCounts = userInput.minBlockCounts;

    // Initialize the first block
    uint64_t blockStart = inputMatches[0].position;
    uint64_t prevPosition = blockStart;
    uint16_t blockCounts = 1;
    uint16_t forwardCount = inputMatches[0].isForward ? 1 : 0;
    uint16_t reverseCount = inputMatches[0].isForward ? 0 : 1;

    auto finalizeBlock = [&](uint64_t endPosition) {
        if (blockCounts >= minBlockCounts) {
            TelomereBlock block;
            block.start = blockStart;
            block.blockLen = (endPosition - blockStart) + patternSize;
            block.blockCounts = blockCounts;

            // p/q assignment
            float forwardRatio = (forwardCount * 100.0f) / blockCounts;
            float reverseRatio = (reverseCount * 100.0f) / blockCounts;
            if (forwardRatio > 60.0f) {
                block.blockLabel = 'p';
            } else if (reverseRatio > 60.0f) {
                block.blockLabel = 'q';
            } else {
                block.blockLabel = 'u';
            }

            telomereBlocks.push_back(block);
        }
    };

    for (size_t i = 1; i < inputMatches.size(); ++i) {
        uint32_t distance = inputMatches[i].position - prevPosition;
        
        if (distance <= mergeDist) {
            // Continue current block
            blockCounts++;
            forwardCount += (inputMatches[i].isForward ? 1 : 0);
            reverseCount += (inputMatches[i].isForward ? 0 : 1);
            prevPosition = inputMatches[i].position;

        } else {
            // Finalize current block
            finalizeBlock(prevPosition);

            // Start new block
            blockStart = inputMatches[i].position;
            prevPosition = blockStart;
            blockCounts = 1;
            forwardCount = (inputMatches[i].isForward ? 1 : 0);
            reverseCount = (inputMatches[i].isForward ? 0 : 1);
        }
    }

    // Finalize the last block
    finalizeBlock(prevPosition);

    return telomereBlocks;
}


std::vector<TelomereBlock> Teloscope::filterTerminalBlocks(const std::vector<TelomereBlock>& blocks) {
    std::vector<TelomereBlock> filteredBlocks;

    // Track 'p' blocks
    TelomereBlock best_p;
    bool has_p1 = false;
    TelomereBlock second_best_p;
    bool has_p2 = false;

    // Track 'q' blocks
    TelomereBlock best_q;
    bool has_q1 = false;
    TelomereBlock second_best_q;
    bool has_q2 = false;

    for (const auto& block : blocks) {
        if (block.blockLen < userInput.minBlockLen || block.blockLabel == 'u') {
            continue; // Exclude spurious blocks
        }

        if (block.blockLabel == 'p') {
            // Update best_p and second_best_p
            if (!has_p1 || block.blockLen > best_p.blockLen) {
                second_best_p = best_p;
                has_p2 = has_p1;
                best_p = block;
                has_p1 = true;
            }
            else if (!has_p2 || block.blockLen > second_best_p.blockLen) {
                second_best_p = block;
                has_p2 = true;
            }
        }
        else if (block.blockLabel == 'q') {
            // Update best_q and second_best_q
            if (!has_q1 || block.blockLen > best_q.blockLen) {
                second_best_q = best_q;
                has_q2 = has_q1;
                best_q = block;
                has_q1 = true;
            }
            else if (!has_q2 || block.blockLen > second_best_q.blockLen) {
                second_best_q = block;
                has_q2 = true;
            }
        }
    }

    // Assign best 'p' and 'q' blocks per path
    if (has_p1 && has_q1) { 
        filteredBlocks.push_back(best_p);
        filteredBlocks.push_back(best_q);
    }
    else { 
        if (has_p1) {
            filteredBlocks.push_back(best_p);
            if (has_p2) {
                filteredBlocks.push_back(second_best_p);
            }
        }
        if (has_q1) {
            filteredBlocks.push_back(best_q);
            if (has_q2) {
                filteredBlocks.push_back(second_best_q);
            }
        }
    }

    // Sort by coords
    std::sort(filteredBlocks.begin(), filteredBlocks.end(),
            [](const TelomereBlock &a, const TelomereBlock &b) {
                return a.start < b.start;
            });

    return filteredBlocks;
}

std::vector<TelomereBlock> Teloscope::filterInterstitialBlocks(
            const std::vector<TelomereBlock>& interstitialBlocks,
            const std::vector<TelomereBlock>& terminalBlocks) {
    // If there are no terminal blocks, everything is interstitial
    if (terminalBlocks.empty()) {
        return interstitialBlocks;
    }

    std::vector<TelomereBlock> filteredBlocks;
    filteredBlocks.reserve(interstitialBlocks.size());

    // Helper lambda to check overlap
    auto overlapsWith = [&](const TelomereBlock &intr, const TelomereBlock &term) {
        uint64_t intrEnd = intr.start + intr.blockLen;
        uint64_t termEnd = term.start + term.blockLen;
        return (intr.start < termEnd) && (intrEnd > term.start);
    };

    // Trim interstitial blocks from the start
    size_t startIdx = 0;
    while (startIdx < interstitialBlocks.size() && overlapsWith(interstitialBlocks[startIdx], terminalBlocks.front())) {
        ++startIdx;
    }

    // Trim interstitial blocks from the end
    size_t endIdx = interstitialBlocks.size();
    while (endIdx > startIdx && overlapsWith(interstitialBlocks[endIdx - 1], terminalBlocks.back())) {
        --endIdx;
    }

    // Collect the remaining interstitial blocks
    for (size_t i = startIdx; i < endIdx; ++i) {
        filteredBlocks.push_back(interstitialBlocks[i]);
    }

    return filteredBlocks;
}

void Teloscope::analyzeWindow(const std::string_view &window, uint32_t windowStart,
                            WindowData& windowData, WindowData& nextOverlapData,
                            SegmentData& segmentData, uint32_t segmentSize, uint32_t absPos) {

    windowData.windowStart = windowStart; // CHECK: Why is this here?
    unsigned short int longestPatternSize = this->trie.getLongestPatternSize();
    uint32_t step = userInput.step;
    uint32_t overlapSize = userInput.windowSize - step;
    uint32_t terminalLimit = userInput.terminalLimit;
    bool computeGC = userInput.outGC;
    bool computeEntropy = userInput.outEntropy;

    // Determine starting index for Trie scanning
    uint32_t startIndex = (windowStart == 0 || overlapSize == 0)
                            ? 0 
                            : std::min(step - longestPatternSize, overlapSize - longestPatternSize);
    int lastCanonicalPos = -1;

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

        for (uint32_t j = i; j < scanLimit; ++j) { // Scan positions until longest pattern
            current = trie.getChild(current, window[j]); // window[j] is a character
            if (!current) break;

            if (current->isEndOfWord) {                
                std::string_view pattern(window.data() + i, (j - i + 1));
                bool isForward = (pattern.size() >= 3 && pattern.compare(0, 3, "CCC") == 0);
                bool isCanonical = (pattern == std::string_view(userInput.canonicalPatterns.first) || 
                                    pattern == std::string_view(userInput.canonicalPatterns.second));
                bool isTerminal = (windowStart + i <= terminalLimit || 
                                    windowStart + i >= segmentSize - terminalLimit);
                float densityGain = static_cast<float>(pattern.size()) / window.size();
                uint32_t matchPos = absPos + windowStart + i; // Keep absolute positions only

                MatchInfo matchInfo;
                matchInfo.position = matchPos;
                matchInfo.isCanonical = isCanonical;
                matchInfo.isForward = isForward;

                // Update windowData
                if (j >= overlapSize || overlapSize == 0 || windowStart == 0) {
                    if (isCanonical) {
                        windowData.canonicalCounts++;
                        windowData.canonicalDensity += densityGain;
                        segmentData.canonicalMatches.push_back(matchInfo);

                        // Check for canonical dimer
                        if (!windowData.hasCanDimer && lastCanonicalPos >= 0 &&
                            (matchPos - lastCanonicalPos) <= userInput.canonicalSize) {
                            windowData.hasCanDimer = true;
                        }
                        lastCanonicalPos = matchPos;  // Update last pos

                    } else {
                        windowData.nonCanonicalCounts++;
                        windowData.nonCanonicalDensity += densityGain;
                    }
                    if (isTerminal) {
                        windowData.winMatches.push_back(matchInfo);
                        if (!isCanonical) {
                            segmentData.nonCanonicalMatches.push_back(matchInfo);
                        }
                    }
                }

                // Update nextOverlapData
                if (i >= step && overlapSize != 0) {
                    if (isCanonical) {
                        nextOverlapData.canonicalCounts++;
                        nextOverlapData.canonicalDensity += densityGain;
                    } else {
                        nextOverlapData.nonCanonicalCounts++;
                        nextOverlapData.nonCanonicalDensity += densityGain;
                    }
                }
            }
        }
    }
}


SegmentData Teloscope::analyzeSegment(std::string &sequence, UserInputTeloscope userInput, uint32_t absPos) {
    uint32_t windowSize = userInput.windowSize;
    uint32_t step = userInput.step;
    uint32_t segmentSize = sequence.size();

    SegmentData segmentData;
    segmentData.segBlocks.reserve(segmentSize / 13 + 1);
    segmentData.canonicalMatches.reserve(segmentSize / 6);
    segmentData.nonCanonicalMatches.reserve(segmentSize / 6);
    segmentData.segMatches.reserve(segmentSize / 6);
    segmentData.windows.reserve((segmentSize - windowSize) / step + 2);

    WindowData prevOverlapData; // Data from previous overlap
    WindowData nextOverlapData; // Data for next overlap

    std::vector<WindowData> windows;
    uint32_t windowStart = 0;
    uint32_t currentWindowSize = std::min(windowSize, segmentSize); // In case first segment is short
    std::string_view windowView(sequence.data(), currentWindowSize);

    while (windowStart < segmentSize) {
        // Prepare and analyze current window
        WindowData windowData = prevOverlapData;

        analyzeWindow(windowView, windowStart, 
                    windowData, nextOverlapData, 
                    segmentData, segmentSize, absPos);

        if (userInput.outGC) { windowData.gcContent = getGCContent(windowData.nucleotideCounts, windowView.size()); }
        if (userInput.outEntropy) { windowData.shannonEntropy = getShannonEntropy(windowData.nucleotideCounts, windowView.size()); }

        // Update windowData
        windowData.windowStart = windowStart + absPos;
        windowData.currentWindowSize = currentWindowSize;
        windows.emplace_back(windowData);

        prevOverlapData = nextOverlapData; // Pass and reset overlap data
        nextOverlapData = WindowData(); // Reset for next iteration

        // Keep all in presence of canonical matches
        if (windowData.hasCanDimer) {
            segmentData.segMatches.insert(segmentData.segMatches.end(),
                                        windowData.winMatches.begin(),
                                        windowData.winMatches.end());
        }

        // Check if we reached the end
        windowStart += step;
        if (windowStart >= segmentSize) {
            break;
        }

        // Prepare next window
        currentWindowSize = std::min(windowSize, segmentSize - windowStart);
        windowView = std::string_view(sequence.data() + windowStart, currentWindowSize);
    }

    // Process "all" and "canonical" matches
    if (segmentData.segMatches.size() >= 2) {
        uint16_t mergeDist = userInput.maxBlockDist;
        segmentData.terminalBlocks = getTelomereBlocks(segmentData.segMatches, mergeDist);
    }

    if (segmentData.canonicalMatches.size() >= 2) {
        uint16_t mergeDist = this->trie.getLongestPatternSize();
        segmentData.interstitialBlocks = getTelomereBlocks(segmentData.canonicalMatches, mergeDist);
    }

    segmentData.windows = windows;
    return segmentData;
}


void Teloscope::writeBEDFile(std::ofstream& windowMetricsFile, std::ofstream& windowRepeatsFile,
                            std::ofstream& canonicalMatchFile, std::ofstream& noncanonicalMatchFile,
                            std::ofstream& terminalBlocksFile, std::ofstream& interstitialBlocksFile) {

    // Write header for window_metrics.tsv
    if (userInput.outEntropy || userInput.outGC) {
        windowMetricsFile << "Header\tStart\tEnd";
        if (userInput.outEntropy) {
            windowMetricsFile << "\tShannonEntropy";
        }
        if (userInput.outGC) {
            windowMetricsFile << "\tGCContent";
        }
        windowMetricsFile << "\n";
    }

    for (const auto& pathData : allPathData) {
        const auto& header = pathData.header;
        const auto& windows = pathData.windows;

        // Write blocks
        for (const auto& block : pathData.terminalBlocks) {
            uint64_t blockEnd = block.start + block.blockLen;
            terminalBlocksFile << header << "\t"
                            << block.start << "\t"
                            << blockEnd << "\t"
                            << block.blockLen << "\t"
                            << block.blockLabel << "\n";
        }

        if (userInput.outITS) {
            for (const auto& block : pathData.interstitialBlocks) {
                uint64_t blockEnd = block.start + block.blockLen;
                interstitialBlocksFile << header << "\t"
                            << block.start << "\t"
                            << blockEnd << "\t"
                            << block.blockLen << "\t"
                            << block.blockLabel << "\n";
            }
        }

        // Write matches to separate files
        if (userInput.outMatches) {
            for (const auto& match : pathData.canonicalMatches) {
                canonicalMatchFile << header << "\t"
                                << match.position << "\t"
                                << (match.position + 6) << "\t"
                                << "canonical" << "\n";
            }

            for (const auto& match : pathData.nonCanonicalMatches) {
                noncanonicalMatchFile << header << "\t"
                                    << match.position << "\t"
                                    << (match.position + 6) << "\t"
                                    << "non-canonical" << "\n";
            }
        }


        // Process window data
        for (const auto& window : windows) {
            totalNWindows++; // Update total window count
            uint32_t windowEnd = window.windowStart + window.currentWindowSize; // Start is already 0-based

            // Write window metrics if enabled
            if (userInput.outEntropy || userInput.outGC) {
                windowMetricsFile << header << "\t" << window.windowStart << "\t" << windowEnd;

                if (userInput.outEntropy) {
                    windowMetricsFile << "\t" << window.shannonEntropy;
                    entropyValues.push_back(window.shannonEntropy); // For summary
                }
                if (userInput.outGC) {
                    windowMetricsFile << "\t" << window.gcContent;
                    gcContentValues.push_back(window.gcContent); // For summary
                }
                windowMetricsFile << "\n";
            }

            // Write repeats data if enabled
            if (userInput.outWinRepeats && window.canonicalCounts > 1) {
                windowRepeatsFile << header << "\t" << window.windowStart << "\t"
                                << windowEnd << "\t"
                                << window.canonicalCounts << "\t"
                                << window.nonCanonicalCounts << "\t"
                                << window.canonicalDensity << "\t"
                                << window.nonCanonicalDensity << "\n";
            }
        }
    }
}


void Teloscope::handleBEDFile() {
    std::ofstream windowMetricsFile;
    std::ofstream windowRepeatsFile;
    std::ofstream canonicalMatchFile;
    std::ofstream noncanonicalMatchFile;
    std::ofstream terminalBlocksFile;
    std::ofstream interstitialBlocksFile;

    std::cout << "Reporting window matches and metrics in BED/BEDgraphs...\n";

    // Open files for writing
    if (userInput.outWinRepeats) {
        windowRepeatsFile.open(userInput.outRoute + "/window_repeats.bedgraph");
    }

    if (userInput.outEntropy || userInput.outGC) {
        windowMetricsFile.open(userInput.outRoute + "/window_metrics.tsv");
    }

    if (userInput.outMatches) {
        canonicalMatchFile.open(userInput.outRoute + "/canonical_matches.bed");
        noncanonicalMatchFile.open(userInput.outRoute + "/noncanonical_matches.bed");
    }

    if (userInput.outITS) {
        interstitialBlocksFile.open(userInput.outRoute + "/interstitial_telomeres.bed");
    }

    terminalBlocksFile.open(userInput.outRoute + "/terminal_telomeres.bed");

    writeBEDFile(windowMetricsFile, windowRepeatsFile, canonicalMatchFile, 
                noncanonicalMatchFile, terminalBlocksFile, interstitialBlocksFile);

    // Close all files once
    if (userInput.outWinRepeats) {
        windowRepeatsFile.close();
    }

    if (userInput.outEntropy || userInput.outGC) {
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
    std::cout << "\n+++Summary Report+++\n";
    std::cout << "Total windows analyzed:\t" << totalNWindows << "\n";

    if (userInput.outEntropy) {
        Stats entropyStats = getStats(entropyValues);
        std::cout << "Max Shannon Entropy:\t" << entropyStats.max << "\n";
        std::cout << "Mean Shannon Entropy:\t" << entropyStats.mean << "\n";
        std::cout << "Median Shannon Entropy:\t" << entropyStats.median << "\n";
        std::cout << "Min Shannon Entropy:\t" << entropyStats.min << "\n";
    }
    
    if (userInput.outGC) {
        Stats gcContentStats = getStats(gcContentValues);
        std::cout << "Max GC Content:\t" << gcContentStats.max << "\n";
        std::cout << "Mean GC Content:\t" << gcContentStats.mean << "\n";
        std::cout << "Median GC Content:\t" << gcContentStats.median << "\n";
        std::cout << "Min GC Content:\t" << gcContentStats.min << "\n";
    }
}