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


std::vector<TelomereBlock> Teloscope::getTelomereBlocks(const std::vector<uint32_t>& inputMatches, uint16_t mergeDist) {
    std::vector<TelomereBlock> telomereBlocks;
    uint16_t patternSize = userInput.canonicalSize;
    uint16_t minBlockCounts = userInput.minBlockCounts;

    if (inputMatches.empty()) {
        return telomereBlocks;
    }

    // Initialize the first block
    uint64_t blockStart = inputMatches[0];
    uint64_t prevPosition = blockStart;
    uint16_t blockCounts = 1;

    auto finalizeBlock = [&](uint64_t endPosition) {
        if (blockCounts >= minBlockCounts) {
            TelomereBlock block;
            block.start = blockStart;
            block.blockLen = (endPosition - blockStart) + patternSize;
            block.blockCounts = blockCounts;
            telomereBlocks.push_back(block);
        }
    };

    for (size_t i = 1; i <= inputMatches.size(); ++i) {
        uint64_t currentPosition;
        uint64_t distance;

        if (i < inputMatches.size()) {
            currentPosition = inputMatches[i];
            distance = currentPosition - prevPosition;
        } else {
            finalizeBlock(prevPosition);
            break;
        }

        if (distance <= mergeDist) {
            blockCounts++; // Extend the block
        } else {
            finalizeBlock(prevPosition);
            blockStart = currentPosition; // Start a new block
            blockCounts = 1;
        }
        prevPosition = currentPosition;
    }

    return telomereBlocks;
}


std::vector<TelomereBlock> Teloscope::filterBlocks(const std::vector<TelomereBlock>& blocks) {
    std::vector<TelomereBlock> filteredBlocks;

    for (const auto& block : blocks) {
        if (block.blockLen >= userInput.minBlockLen) { // Filter by length
            filteredBlocks.push_back(block);
        }
    }

    if (filteredBlocks.size() > 2) {  // Keep the two largest blocks
        std::nth_element(filteredBlocks.begin(), filteredBlocks.begin() + 2, filteredBlocks.end(),
            [](const TelomereBlock& a, const TelomereBlock& b) {
                return a.blockLen > b.blockLen;
            });
        filteredBlocks.resize(2);
    }

    // Ensure the blocks are ordered by their start positions
    std::sort(filteredBlocks.begin(), filteredBlocks.end(), [](const TelomereBlock& a, const TelomereBlock& b) {
        return a.start < b.start;
    });

    return filteredBlocks;
}


void Teloscope::analyzeWindow(const std::string &window, uint32_t windowStart,
                            WindowData& windowData, WindowData& nextOverlapData,
                            SegmentData& segmentData, uint32_t segmentSize) {

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
            if (i >= overlapSize || overlapSize == 0 || windowStart == 0) {
                switch (window[i]) {
                    case 'A': windowData.nucleotideCounts[0]++; break;
                    case 'C': windowData.nucleotideCounts[1]++; break;
                    case 'G': windowData.nucleotideCounts[2]++; break;
                    case 'T': windowData.nucleotideCounts[3]++; break;
                }
            }
            if (i >= step && overlapSize != 0) {
                switch (window[i]) {
                    case 'A': nextOverlapData.nucleotideCounts[0]++; break;
                    case 'C': nextOverlapData.nucleotideCounts[1]++; break;
                    case 'G': nextOverlapData.nucleotideCounts[2]++; break;
                    case 'T': nextOverlapData.nucleotideCounts[3]++; break;
                }
            }
        }

        // Pattern matching using Trie
        auto current = trie.getRoot();
        uint32_t scanLimit = std::min(i + longestPatternSize, static_cast<uint32_t>(window.size()));

        for (uint32_t j = i; j < scanLimit; ++j) { // Scan positions until longest pattern
            current = trie.getChild(current, window[j]); // window[j] is a character
            if (!current) break;

            if (current->isEndOfWord) {
                std::string pattern = window.substr(i, j - i + 1);
                bool isCanonical = (pattern == userInput.canonicalPatterns.first || 
                                    pattern == userInput.canonicalPatterns.second); // Check canonical patterns
                bool isTerminal = (windowStart + i <= terminalLimit || 
                                    windowStart + i >= segmentSize - terminalLimit);
                float densityGain = static_cast<float>(pattern.size()) / window.size();
                uint32_t matchPos = windowStart + i;

                // Update windowData
                if (j >= overlapSize || overlapSize == 0 || windowStart == 0) {
                    if (isCanonical) {
                        windowData.canonicalCounts++;
                        windowData.canonicalDensity += densityGain;
                        segmentData.canonicalMatches.push_back(matchPos);

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
                        windowData.winMatches.push_back(matchPos);
                        if (!isCanonical) {
                            segmentData.nonCanonicalMatches.push_back(matchPos);
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


SegmentData Teloscope::analyzeSegment(std::string &sequence, UserInputTeloscope userInput, uint64_t absPos) {
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
    std::unordered_map<std::string, std::vector<TelomereBlock>> mergedBlocks;
    uint32_t windowStart = 0;
    uint32_t currentWindowSize = std::min(windowSize, segmentSize); // In case first segment is short
    std::string window = sequence.substr(0, currentWindowSize);

    while (windowStart < segmentSize) {
        // Prepare and analyze current window
        WindowData windowData = prevOverlapData;
        analyzeWindow(window, windowStart, windowData, nextOverlapData, segmentData, segmentSize);

        if (userInput.outGC) { windowData.gcContent = getGCContent(windowData.nucleotideCounts, window.size()); }
        if (userInput.outEntropy) { windowData.shannonEntropy = getShannonEntropy(windowData.nucleotideCounts, window.size()); }

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

        // Prepare next window
        windowStart += step;
        if (windowStart >= segmentSize) {
            break;
        }

        // Recycle the overlapping string sequence
        currentWindowSize = std::min(windowSize, segmentSize - windowStart);
        if (currentWindowSize == windowSize) {
            window = window.substr(step) + sequence.substr(windowStart + windowSize - step, step);
        } else {
            window = sequence.substr(windowStart, currentWindowSize); // Last window has a shorter size
        }
    }

    // Process "all" and "canonical" matches
    if (segmentData.segMatches.size() >= 2) {
        uint16_t mergeDist = userInput.maxBlockDist;
        auto segBlocks = getTelomereBlocks(segmentData.segMatches, mergeDist);
        segmentData.mergedBlocks["all"] = segBlocks;
    }

    if (segmentData.canonicalMatches.size() >= 2) {
        uint16_t mergeDist = this->trie.getLongestPatternSize();
        auto canonicalBlocks = getTelomereBlocks(segmentData.canonicalMatches, mergeDist);
        segmentData.mergedBlocks["canonical"] = canonicalBlocks;
    }

    segmentData.windows = windows;
    return segmentData;
}


void Teloscope::writeBEDFile(std::ofstream& windowMetricsFile, std::ofstream& windowRepeatsFile,
                            std::ofstream& canonicalMatchFile, std::ofstream& noncanonicalMatchFile,
                            std::ofstream& allBlocksFile, std::ofstream& canonicalBlocksFile) {

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

        // Write telomere blocks
        for (const auto& [groupName, blocks] : pathData.mergedBlocks) {
            std::ofstream* outputFile = nullptr;

            // By group
            if (groupName == "all") {
                outputFile = &allBlocksFile;
            } else if (groupName == "canonical") {
                outputFile = &canonicalBlocksFile;
            }

            if (outputFile) {
                for (const auto& block : blocks) {
                    uint64_t blockEnd = block.start + block.blockLen;
                    *outputFile << header << "\t" << block.start << "\t" << blockEnd << "\t" << block.blockLen << "\n";
                }
            }
        }

        // Write matches to separate files.
        for (const auto& pos : pathData.canonicalMatches) {
            canonicalMatchFile << header << "\t"
                                << pos << "\t"
                                << pos + 6 << "\t"
                                << "canonical" << "\n";
        }

        for (const auto& pos : pathData.nonCanonicalMatches) {
            noncanonicalMatchFile << header << "\t"
                                    << pos << "\t"
                                    << pos + 6 << "\t"
                                    << "non-canonical" << "\n";
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
            if (window.canonicalCounts > 1) {
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
    std::ofstream allBlocksFile;
    std::ofstream canonicalBlocksFile;

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



    allBlocksFile.open(userInput.outRoute + "/telomere_blocks_all.bed");
    canonicalBlocksFile.open(userInput.outRoute + "/telomere_blocks_canonical.bed");

    writeBEDFile(windowMetricsFile, windowRepeatsFile, canonicalMatchFile, 
                noncanonicalMatchFile, allBlocksFile, canonicalBlocksFile);

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



    allBlocksFile.close();
    canonicalBlocksFile.close();
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