#include <iostream>
#include <fstream>
#include <sstream> // check
#include <stdint.h> // what's this for?
#include <vector>
#include <algorithm> // removeCarriageReturns
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


float Teloscope::getMean(const std::vector<float>& values) {
    if (values.empty()) return 0.0;
    float sum = std::accumulate(values.begin(), values.end(), 0.0);
    return sum / values.size();
}


float Teloscope::getMedian(std::vector<float> values) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    size_t size = values.size();
    if (size % 2 == 0) {
        return (values[size / 2 - 1] + values[size / 2]) / 2;
    } else {
        return values[size / 2];
    }
}


float Teloscope::getMin(const std::vector<float> values) {
    if (values.empty()) return 0.0;
    return *std::min_element(values.begin(), values.end());
}


float Teloscope::getMax(const std::vector<float> values) {
    if (values.empty()) return 0.0;
    return *std::max_element(values.begin(), values.end());
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



void Teloscope::analyzeWindow(const std::string &window, uint32_t windowStart,
                            WindowData& windowData, WindowData& nextOverlapData,
                            SegmentData& segmentData, uint32_t segmentSize) {
// void Teloscope::analyzeWindow(const std::string &window, uint32_t windowStart, WindowData& windowData, WindowData& nextOverlapData, SegmentData& segmentData) {

    windowData.windowStart = windowStart; // CHECK: Why is this here?
    unsigned short int longestPatternSize = this->trie.getLongestPatternSize();
    uint32_t overlapSize = userInput.windowSize - userInput.step;
    uint32_t terminalLimit = userInput.terminalLimit;

    // Determine starting index for Trie scanning
    uint32_t startIndex = (windowStart == 0 || overlapSize == 0) ? 0 : std::min(userInput.step - longestPatternSize, overlapSize - longestPatternSize);

    for (uint32_t i = startIndex; i < window.size(); ++i) { // Nucleotide iterations

        if (userInput.modeGC || userInput.modeEntropy) {
            if (i >= overlapSize || overlapSize == 0 || windowStart == 0) {
                switch (window[i]) {
                    case 'A': windowData.nucleotideCounts[0]++; break;
                    case 'C': windowData.nucleotideCounts[1]++; break;
                    case 'G': windowData.nucleotideCounts[2]++; break;
                    case 'T': windowData.nucleotideCounts[3]++; break;
                }
            }
            if (i >= userInput.step && userInput.windowSize != userInput.step) {
                switch (window[i]) {
                    case 'A': nextOverlapData.nucleotideCounts[0]++; break;
                    case 'C': nextOverlapData.nucleotideCounts[1]++; break;
                    case 'G': nextOverlapData.nucleotideCounts[2]++; break;
                    case 'T': nextOverlapData.nucleotideCounts[3]++; break;
                }
            }
        }

        if (userInput.modeMatch) { // Pattern matching using Trie
            auto current = trie.getRoot();
            uint32_t scanLimit = std::min(i + longestPatternSize, static_cast<uint32_t>(window.size()));

            for (uint32_t j = i; j < scanLimit; ++j) { // Scan positions until longest pattern

                if (!trie.hasChild(current, window[j])) break;  
                current = trie.getChild(current, window[j]); // window[j] is a character

                if (current->isEndOfWord) {
                    std::string pattern = window.substr(i, j - i + 1);
                    bool isCanonical = (pattern == userInput.canonicalPatterns.first || pattern == userInput.canonicalPatterns.second); // Check canonical patterns
                    bool isTerminal = (windowStart + i <= terminalLimit || windowStart + i >= segmentSize - terminalLimit);
                    uint8_t patternSize = pattern.size();

                    // Update windowData from prevOverlapData
                    if (j >= overlapSize || overlapSize == 0 || windowStart == 0) {
                        float densityGain = static_cast<float>(patternSize) / window.size();
                        uint32_t matchPos = windowStart + i;

                        if (isCanonical) {
                            windowData.canonicalCounts++;
                            windowData.canonicalDensity += densityGain;
                            segmentData.canonicalMatches.push_back(matchPos);
                        } else {
                            windowData.nonCanonicalCounts++;
                            windowData.nonCanonicalDensity += densityGain;
                        }
                        if (isTerminal) {
                            segmentData.segMatches.push_back(matchPos);
                            if (!isCanonical) {
                                segmentData.nonCanonicalMatches.push_back(matchPos);
                            }
                        }
                        // windowData.hDistances.push_back(userInput.hammingDistances[pattern]);
                        // windowData.winHDistance += userInput.hammingDistances[pattern];
                    }

                    // Update nextOverlapData
                    if (i >= userInput.step && overlapSize != 0 ) {
                        if (isCanonical) {
                            nextOverlapData.canonicalCounts++;
                            nextOverlapData.canonicalDensity += static_cast<float>(patternSize) / window.size();
                        } else {
                            nextOverlapData.nonCanonicalCounts++;
                            nextOverlapData.nonCanonicalDensity += static_cast<float>(patternSize) / window.size();
                        }
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
        // analyzeWindow(window, windowStart, windowData, nextOverlapData, segmentData);
        analyzeWindow(window, windowStart, windowData, nextOverlapData, segmentData, segmentSize);

        if (userInput.modeGC) { windowData.gcContent = getGCContent(windowData.nucleotideCounts, window.size()); }
        if (userInput.modeEntropy) { windowData.shannonEntropy = getShannonEntropy(windowData.nucleotideCounts, window.size()); }

        // Update windowData
        windowData.windowStart = windowStart + absPos;
        windowData.currentWindowSize = currentWindowSize;
        if (userInput.keepWindowData) { windows.emplace_back(windowData); } // User defines if windowData is stored

        // Pass and reset overlap data
        prevOverlapData = nextOverlapData;
        nextOverlapData = WindowData(); // Reset for next iteration

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
    if (userInput.keepWindowData && (userInput.modeEntropy || userInput.modeGC)) {
        windowMetricsFile << "Header\tStart\tEnd";
        if (userInput.modeEntropy) {
            windowMetricsFile << "\tShannonEntropy";
        }
        if (userInput.modeGC) {
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
                    *outputFile << header << "\t" << block.start << "\t" << blockEnd << "\t" << groupName << "\n";
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

        if (!userInput.keepWindowData) { // If windowData is not stored, return
            continue;
        }

        // Process window data
        for (const auto& window : windows) {
            totalNWindows++; // Update total window count
            uint32_t windowEnd = window.windowStart + window.currentWindowSize; // Start is already 0-based

            // Write window metrics if enabled
            if (userInput.modeEntropy || userInput.modeGC) {
                windowMetricsFile << header << "\t" << window.windowStart << "\t" << windowEnd;

                if (userInput.modeEntropy) {
                    windowMetricsFile << "\t" << window.shannonEntropy;
                    entropyValues.push_back(window.shannonEntropy); // For summary
                }
                if (userInput.modeGC) {
                    windowMetricsFile << "\t" << window.gcContent;
                    gcContentValues.push_back(window.gcContent); // For summary
                }
                windowMetricsFile << "\n";
            }

            // Write repeats data if enabled
            if (userInput.modeMatch && window.canonicalCounts > 1) {
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
    if (userInput.keepWindowData && (userInput.modeEntropy || userInput.modeGC)) {
        windowMetricsFile.open(userInput.outRoute + "/window_metrics.tsv");
    }

    if (userInput.keepWindowData && userInput.modeMatch) {
        canonicalMatchFile.open(userInput.outRoute + "/canonical_matches.bed");
        noncanonicalMatchFile.open(userInput.outRoute + "/noncanonical_matches.bed");
        windowRepeatsFile.open(userInput.outRoute + "/window_repeats.bedgraph");
    }

    allBlocksFile.open(userInput.outRoute + "/telomere_blocks_all.bed");
    canonicalBlocksFile.open(userInput.outRoute + "/telomere_blocks_canonical.bed");

    writeBEDFile(windowMetricsFile, windowRepeatsFile, canonicalMatchFile, 
                noncanonicalMatchFile, allBlocksFile, canonicalBlocksFile);

    // Close all files once
    if (userInput.modeEntropy || userInput.modeGC) {
        windowMetricsFile.close();
    }
    if (userInput.modeMatch) {
        canonicalMatchFile.close();
        noncanonicalMatchFile.close();
        windowRepeatsFile.close();
    }

    allBlocksFile.close();
    canonicalBlocksFile.close();
}


void Teloscope::printSummary() {
    if (!userInput.keepWindowData) { // If windowData is not stored, skip
        return;
    }

    std::cout << "\n+++Summary Report+++\n";
    std::cout << "Total windows analyzed:\t" << totalNWindows << "\n";
    // Print the total canonical and non-canonical matches per path - PENDING
    // For each pattern, print the path header with the highest number of matches - PENDING
    // For each pattern, print the path header with the lowest number of matches - PENDING
    std::cout << "Max Shannon Entropy:\t" << getMax(entropyValues) << "\n";
    std::cout << "Mean Shannon Entropy:\t" << getMean(entropyValues) << "\n";
    std::cout << "Median Shannon Entropy:\t" << getMedian(entropyValues) << "\n";
    std::cout << "Min Shannon Entropy:\t" << getMin(entropyValues) << "\n";

    std::cout << "Max GC Content:\t" << getMax(gcContentValues) << "\n";
    std::cout << "Mean GC Content:\t" << getMean(gcContentValues) << "\n";
    std::cout << "Median GC Content:\t" << getMedian(gcContentValues) << "\n";
    std::cout << "Min GC Content:\t" << getMin(gcContentValues) << "\n";
}