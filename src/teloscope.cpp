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


float Teloscope::getShannonEntropy(const std::unordered_map<char, uint32_t>& nucleotideCounts, uint32_t windowSize) {
    float entropy = 0.0;
    for (auto &[nucleotide, count] : nucleotideCounts) {
        if (count > 0) {
            float probability = static_cast<float>(count) / windowSize;
            entropy -= probability * std::log2(probability);
        }
    }
    return entropy;
}


float Teloscope::getGCContent(const std::unordered_map<char, uint32_t>& nucleotideCounts, uint32_t windowSize) {
    uint32_t gcCount = nucleotideCounts.at('G') + nucleotideCounts.at('C');
    return float(gcCount) / windowSize * 100.0;
}


void Teloscope::getPatternDensities(WindowData& windowData, uint32_t windowSize) {
    for (auto &entry : windowData.patternMap) {
        auto &pattern = entry.first;
        auto &data = entry.second;
        data.density = static_cast<float>(data.count * pattern.size()) / windowSize;
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


std::vector<TelomereBlock> Teloscope::getTelomereBlocks(const std::vector<uint32_t>& inputMatches, uint64_t windowStart) {
    std::vector<TelomereBlock> winBlocks;
    uint16_t patternSize = 6;
    uint16_t D = this->trie.getLongestPatternSize(); // D is set to longestPatternSize

    if (inputMatches.empty()) {
        return winBlocks; // No matches to process
    }

    // Initialize the first block
    uint64_t blockStart = windowStart + inputMatches[0];
    uint64_t prevPosition = blockStart;
    uint16_t blockCounts = 1;
    
    // Helper function
    auto finalizeBlock = [&](uint64_t endPosition) {
        if (blockCounts >= 2) {
            TelomereBlock block;
            block.start = blockStart;
            block.blockLen = (endPosition - blockStart) + patternSize;
            winBlocks.push_back(block);
        }
    };

    for (size_t i = 1; i <= inputMatches.size(); ++i) {
        uint64_t currentPosition;
        uint64_t distance;

        if (i < inputMatches.size()) {
            currentPosition = windowStart + inputMatches[i];
            distance = currentPosition - prevPosition;
        } else {
            currentPosition = 0;
            distance = D + 1;    // Finalize last block
        }

        if (distance <= D) {
            blockCounts++; // Extend the block
        } else {
            finalizeBlock(prevPosition); // Finalize the current block

            if (i < inputMatches.size()) {
                blockStart = currentPosition; // Start a new block
                blockCounts = 1;
            }
        }
        prevPosition = currentPosition;
    }

    return winBlocks;
}


std::vector<TelomereBlock> Teloscope::mergeTelomereBlocks(const std::vector<TelomereBlock>& winBlocks) {
    std::vector<TelomereBlock> mergedBlocks;

    if (winBlocks.empty()) {
        return mergedBlocks; // No blocks to merge
    }

    // Initialize the first block as the current merged block
    TelomereBlock currentBlock = winBlocks[0];
    uint16_t D = this->trie.getLongestPatternSize(); // Use D as the merging distance threshold

    for (size_t i = 1; i < winBlocks.size(); ++i) {
        const TelomereBlock& nextBlock = winBlocks[i];
        uint64_t currentEnd = currentBlock.start + currentBlock.blockLen;
        uint64_t distance = nextBlock.start - currentEnd;

        if (distance <= 10 * D || distance > 5) {
            uint64_t newEnd = nextBlock.start + nextBlock.blockLen;
            currentBlock.blockLen = newEnd - currentBlock.start;
            
        } else {
            mergedBlocks.push_back(currentBlock);
            currentBlock = nextBlock; // Start a new current block
        }
    }

    // Add the last current block to merged blocks
    mergedBlocks.push_back(currentBlock);

    return mergedBlocks;
}


void Teloscope::analyzeWindow(const std::string &window, uint32_t windowStart, WindowData& windowData, WindowData& nextOverlapData) {
    windowData.windowStart = windowStart; // CHECK: Why is this here?
    unsigned short int longestPatternSize = this->trie.getLongestPatternSize();
    uint32_t overlapSize = userInput.windowSize - userInput.step;

    // Determine starting index for Trie scanning
    uint32_t startIndex = (windowStart == 0 || overlapSize == 0) ? 0 : std::min(userInput.step - longestPatternSize, overlapSize - longestPatternSize);

    for (uint32_t i = startIndex; i < window.size(); ++i) { // Nucleotide iterations

        if (userInput.modeGC || userInput.modeEntropy) {
            if (i >= overlapSize || overlapSize == 0 || windowStart == 0) {
                windowData.nucleotideCounts[window[i]]++; // For whole window
            }
            if (i >= userInput.step && userInput.windowSize != userInput.step) {
                nextOverlapData.nucleotideCounts[window[i]]++; // For next overlap
            }
        }

        if (userInput.modeMatch) { // Pattern matching using Trie
            auto current = trie.getRoot();
            uint32_t scanLimit = std::min(i + longestPatternSize, static_cast<uint32_t>(window.size()));

            for (uint32_t j = i; j < scanLimit; ++j) { // Only scan positions in range of patterns

                if (!trie.hasChild(current, window[j])) break;  
                current = trie.getChild(current, window[j]); // window[j] is a character

                if (current->isEndOfWord) {
                    std::string pattern = window.substr(i, j - i + 1);
                    bool isCanonical = (pattern == userInput.canonicalPatterns.first || pattern == userInput.canonicalPatterns.second); // Check canonical patterns

                    // Update windowData from prevOverlapData
                    if (j >= overlapSize || overlapSize == 0 || windowStart == 0 ) {
                        windowData.patternMap[pattern].count++;
                        isCanonical ? windowData.canonicalCounts++ : windowData.nonCanonicalCounts++;
                        windowData.windowCounts++;

                        windowData.patternMap[pattern].patMatches.push_back(i);
                        isCanonical ? windowData.canonicalMatches.push_back(i) : windowData.nonCanonicalMatches.push_back(i);
                        windowData.windowMatches.push_back(i); // Ordered by design
                    }

                    // Update nextOverlapData
                    if (i >= userInput.step && overlapSize != 0 ) {
                        nextOverlapData.patternMap[pattern].count++;
                        isCanonical ? nextOverlapData.canonicalCounts++ : nextOverlapData.nonCanonicalCounts++;
                    }
                }
            }
        }
    }
}


SegmentData Teloscope::analyzeSegment(std::string &sequence, UserInputTeloscope userInput, uint64_t absPos) {
    uint32_t windowSize = userInput.windowSize;
    uint32_t step = userInput.step;
    uint32_t sequenceSize = sequence.size();

    WindowData prevOverlapData; // Data from previous overlap
    WindowData nextOverlapData; // Data for next overlap

    std::vector<WindowData> windows;
    uint32_t windowStart = 0;
    uint32_t currentWindowSize = std::min(windowSize, static_cast<uint32_t>(sequenceSize)); // In case first segment is short
    std::string window = sequence.substr(0, currentWindowSize);

    std::unordered_map<std::string, std::vector<TelomereBlock>> segmentBlocks = {
        {"all", {}},
        {"canonical", {}},
        {"non-canonical", {}}
    };
    std::unordered_map<std::string, std::vector<TelomereBlock>> mergedBlocks = {
        {"all", {}},
        {"canonical", {}},
        {"non-canonical", {}}
    };

    while (windowStart < sequenceSize) {

        // Prepare and analyze current window
        WindowData windowData = prevOverlapData;
        analyzeWindow(window, windowStart, windowData, nextOverlapData);

        if (userInput.modeGC) {windowData.gcContent = getGCContent(windowData.nucleotideCounts, window.size());}
        if (userInput.modeEntropy) {windowData.shannonEntropy = getShannonEntropy(windowData.nucleotideCounts, window.size());}
        if (userInput.modeMatch) {getPatternDensities(windowData, window.size());}

        // Update windowData
        windowData.windowStart = windowStart + absPos;
        windowData.currentWindowSize = currentWindowSize;
        if (userInput.keepWindowData) {windows.emplace_back(windowData);} // User defines if windowData is stored

        // Collect telomere blocks for each pattern group
        std::unordered_map<std::string, std::vector<uint32_t>> patternMatches = {
            {"all", windowData.windowMatches},
            {"canonical", windowData.canonicalMatches},
            {"non-canonical", windowData.nonCanonicalMatches}
        };

        for (const auto& [groupName, matches] : patternMatches) {
            if (windowData.canonicalCounts >= 2 || windowData.nonCanonicalCounts >= 4) { // JACK: Add to user cutoffs
                auto winBlocks = getTelomereBlocks(matches, windowData.windowStart);
                segmentBlocks[groupName].insert(segmentBlocks[groupName].end(), winBlocks.begin(), winBlocks.end());
            }
        }

        // Pass and reset overlap data
        prevOverlapData = nextOverlapData;
        nextOverlapData = WindowData(); // Reset for next iteration

        // Prepare next window
        windowStart += step;
        if (windowStart >= sequenceSize) {
            break;
        }

        // Recycle the overlapping string sequence
        currentWindowSize = std::min(windowSize, static_cast<uint32_t>(sequenceSize - windowStart)); // CHECK

        if (currentWindowSize == windowSize) {
            window = window.substr(step) + sequence.substr(windowStart + windowSize - step, step);
        } else {
            window = sequence.substr(windowStart, currentWindowSize); // Last window has a shorter size
        }
    }

    // After processing all windows, merge telomere blocks for each pattern group
    for (const auto& [groupName, blocks] : segmentBlocks) {
        if (!blocks.empty()) {
            mergedBlocks[groupName] = mergeTelomereBlocks(blocks);
        }
    }

    SegmentData segmentData;
    segmentData.windows = windows;
    segmentData.mergedBlocks = mergedBlocks;

    return segmentData;

}


// void Teloscope::writeBEDFile(std::ofstream& shannonFile, std::ofstream& gcContentFile,
//                 std::unordered_map<std::string, std::ofstream>& patternMatchFiles,
//                 std::unordered_map<std::string, std::ofstream>& patternCountFiles,
//                 std::unordered_map<std::string, std::ofstream>& patternDensityFiles,
//                 std::ofstream& telomereBlocksFile) {

    // for (const auto& pathData : allPathData) {
    //     // const auto& seqPos = pathData.seqPos;
    //     const auto& header = pathData.header;
    //     const auto& windows = pathData.windows;

    //     // Process telomere blocks
    //     for (const auto& [groupName, blocks] : pathData.mergedBlocks) {
    //         for (const auto& block : blocks) {
    //             uint64_t blockEnd = block.start + block.blockLen;
    //             telomereBlocksFile << header << "\t" << block.start << "\t" << blockEnd << "\t" << groupName << "\n";
    //         }
    //     }

void Teloscope::writeBEDFile(std::ofstream& shannonFile, std::ofstream& gcContentFile,
                            std::unordered_map<std::string, std::ofstream>& patternMatchFiles,
                            std::unordered_map<std::string, std::ofstream>& patternCountFiles,
                            std::unordered_map<std::string, std::ofstream>& patternDensityFiles,
                            std::ofstream& telomereBlocksAllFile,
                            std::ofstream& telomereBlocksCanonicalFile,
                            std::ofstream& telomereBlocksNonCanonicalFile) {


    for (const auto& pathData : allPathData) {
        const auto& header = pathData.header;
        const auto& windows = pathData.windows;

        // Process telomere blocks
        for (const auto& [groupName, blocks] : pathData.mergedBlocks) {
            std::ofstream* outputFile = nullptr;

            // By group
            if (groupName == "all") {
                outputFile = &telomereBlocksAllFile;
            } else if (groupName == "canonical") {
                outputFile = &telomereBlocksCanonicalFile;
            } else if (groupName == "non-canonical") {
                outputFile = &telomereBlocksNonCanonicalFile;
            }

            if (outputFile) {
                for (const auto& block : blocks) {
                    uint64_t blockEnd = block.start + block.blockLen;
                    *outputFile << header << "\t" << block.start << "\t" << blockEnd << "\t" << groupName << "\n";
                }
            }
        }

        if (!userInput.keepWindowData) { // If windowData is not stored, return
            continue;
        }

        // Process window data
        for (const auto& window : windows) {
            totalNWindows++; // Update total window count
            uint32_t windowEnd = window.windowStart + window.currentWindowSize; // Start is already 0-based

            // Write window Shannon entropy if enabled
            if (userInput.modeEntropy) {
                shannonFile << header << "\t" << window.windowStart << "\t"
                            << windowEnd << "\t"
                            << window.shannonEntropy << "\n";
                entropyValues.push_back(window.shannonEntropy); // Update entropy values
            }

            // Write window GC content if enabled
            if (userInput.modeGC) {
                gcContentFile << header << "\t" << window.windowStart << "\t"
                            << windowEnd << "\t"
                            << window.gcContent << "\n";
                gcContentValues.push_back(window.gcContent);
            }

            // Write pattern data if enabled
            if (userInput.modeMatch) {
                for (const auto& [pattern, data] : window.patternMap) {
                    for (auto pos : data.patMatches) {
                        patternMatchFiles[pattern] << header << "\t"
                                                << window.windowStart + pos << "\t"
                                                << window.windowStart + pos + pattern.length() << "\t" // Start is already 0-based
                                                << pattern << "\n";
                        patternCounts[pattern]++; // Update total pattern counts
                    }

                    patternCountFiles[pattern] << header << "\t" << window.windowStart << "\t"
                                            << windowEnd << "\t"
                                            << data.count << "\n";
                                            
                    patternDensityFiles[pattern] << header << "\t" << window.windowStart << "\t"
                                                << windowEnd << "\t"
                                                << data.density << "\n";
                }
            }
        }
    }
}

void Teloscope::handleBEDFile() {
    std::ofstream shannonFile;
    std::ofstream gcContentFile;
    std::unordered_map<std::string, std::ofstream> patternMatchFiles; // Jack: replace with vector to reduce cache locality?
    std::unordered_map<std::string, std::ofstream> patternCountFiles;
    std::unordered_map<std::string, std::ofstream> patternDensityFiles;
    std::ofstream telomereBlocksAllFile;
    std::ofstream telomereBlocksCanonicalFile;
    std::ofstream telomereBlocksNonCanonicalFile;
    std::cout << "Reporting window matches and metrics in BED/BEDgraphs...\n";

    // Open files for writing
    if (userInput.modeEntropy) {
        shannonFile.open(userInput.outRoute + "/shannonEntropy.bedgraph");
    }

    if (userInput.modeGC) {
        gcContentFile.open(userInput.outRoute + "/gcContent.bedgraph");
    }

    if (userInput.modeMatch) {

        for (const auto& pattern : userInput.patterns) {
            patternMatchFiles[pattern].open(userInput.outRoute + "/" + pattern + "_matches.bed");
            patternCountFiles[pattern].open(userInput.outRoute + "/" + pattern + "_count.bedgraph");
            patternDensityFiles[pattern].open(userInput.outRoute + "/" + pattern + "_density.bedgraph");
        }
    }

    telomereBlocksAllFile.open(userInput.outRoute + "/telomere_blocks_all.bed");
    telomereBlocksCanonicalFile.open(userInput.outRoute + "/telomere_blocks_canonical.bed");
    telomereBlocksNonCanonicalFile.open(userInput.outRoute + "/telomere_blocks_noncanonical.bed");

    // Pass the files to writeBEDFile
    writeBEDFile(shannonFile, gcContentFile, patternMatchFiles, patternCountFiles, patternDensityFiles,
                telomereBlocksAllFile, telomereBlocksCanonicalFile, telomereBlocksNonCanonicalFile);

    // Close all files once
    if (userInput.modeEntropy) {
        shannonFile.close();
    }
    if (userInput.modeGC) {
        gcContentFile.close();
    }
    if (userInput.modeMatch) {
        for (auto& [pattern, file] : patternMatchFiles) {
            file.close();
        }
        for (auto& [pattern, file] : patternCountFiles) {
            file.close();
        }
        for (auto& [pattern, file] : patternDensityFiles) {
            file.close();
        }
    }

    telomereBlocksAllFile.close();
    telomereBlocksCanonicalFile.close();
    telomereBlocksNonCanonicalFile.close();
}

void Teloscope::printSummary() {
    std::cout << "\n+++Summary Report+++\n";
    std::cout << "Total windows analyzed:\t" << totalNWindows << "\n";
    std::cout << "Total input patterns found:\n";
    for (const auto& [pattern, count] : patternCounts) {
        std::cout << "Pattern:\t" << pattern << "\t" << count << "\n";
    }

    // For each pattern, print the path header with the highest number of matches - PENDING
    // For each pattern, print the path header with the lowest number of matches - PENDING
    if (userInput.keepWindowData) {
        std::cout << "Max Shannon Entropy:\t" << getMax(entropyValues) << "\n";
        std::cout << "Mean Shannon Entropy:\t" << getMean(entropyValues) << "\n";
        std::cout << "Median Shannon Entropy:\t" << getMedian(entropyValues) << "\n";
        std::cout << "Min Shannon Entropy:\t" << getMin(entropyValues) << "\n";

        std::cout << "Max GC Content:\t" << getMax(gcContentValues) << "\n";
        std::cout << "Mean GC Content:\t" << getMean(gcContentValues) << "\n";
        std::cout << "Median GC Content:\t" << getMedian(gcContentValues) << "\n";
        std::cout << "Min GC Content:\t" << getMin(gcContentValues) << "\n";
    }
}