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


void Teloscope::insertWindowData(unsigned int seqPos, const std::string& header, std::vector<WindowData>& pathWindows) {
    if (userInput.keepWindowData) {
        allWindows.push_back(std::make_tuple(seqPos, header, pathWindows));
    }
}


void Teloscope::sortWindowsBySeqPos() {
    std::sort(allWindows.begin(), allWindows.end(), [](const auto& one, const auto& two) {
        return std::get<0>(one) < std::get<0>(two);
    });
}


void Teloscope::analyzeWindow(const std::string &window, uint32_t windowStart, WindowData& windowData, WindowData& nextOverlapData) {
    windowData.windowStart = windowStart; // CHECK: Why is this here?
    unsigned short int longestPatternSize = this->trie.getLongestPatternSize();
    uint32_t overlapSize = userInput.windowSize - userInput.step; // Overlap starts at this index

    // Determine starting index for Trie scan
    uint32_t startIndex = 0;
    if (windowStart == 0) {
        startIndex = 0;
    } else if (userInput.step < overlapSize) {
        startIndex = userInput.step - longestPatternSize; // To capture patterns missed in the previous window
    } else {
        startIndex = overlapSize - longestPatternSize; // To capture patterns missed in the previous window
    }
    // uint32_t startIndex = (windowStart == 0) ? 0 : std::min(userInput.step, overlapSize); // CHECK: Test whether this is faster

    for (uint32_t i = startIndex; i < window.size(); ++i) { // Nucleotide iterations

        if (userInput.modeGC || userInput.modeEntropy) {
            if (i >= overlapSize || userInput.windowSize == userInput.step || windowStart == 0) {
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

                    // Update windowData from prevOverlapData
                    if (j >= overlapSize || userInput.windowSize == userInput.step || windowStart == 0 ) {
                        windowData.patternMap[pattern].count++;
                        windowData.patternMap[pattern].wMatches.push_back(i);
                    }

                    // Update nextOverlapData from steps
                    if (i >= userInput.step && userInput.windowSize != userInput.step ) {
                        nextOverlapData.patternMap[pattern].count++;
                    }
                }
            }
        }
    }

    if (userInput.modeGC) {
        windowData.gcContent = getGCContent(windowData.nucleotideCounts, window.size());
    }

    if (userInput.modeEntropy) {
        windowData.shannonEntropy = getShannonEntropy(windowData.nucleotideCounts, window.size());
    }

    if (userInput.modeMatch) {
        getPatternDensities(windowData, window.size());
    }
}


std::vector<WindowData> Teloscope::analyzeSegment(std::string &sequence, UserInputTeloscope userInput, uint64_t absPos) {
    uint32_t windowSize = userInput.windowSize;
    uint32_t step = userInput.step;
    uint32_t sequenceSize = sequence.size();

    WindowData prevOverlapData; // Data from previous overlap
    WindowData nextOverlapData; // Data for next overlap

    std::vector<WindowData> windows;
    uint32_t windowStart = 0;
    uint32_t currentWindowSize = std::min(userInput.windowSize, static_cast<uint32_t>(sequence.size())); // In case first segment is short
    std::string window = sequence.substr(0, currentWindowSize);

    while (windowStart < sequenceSize) {

        // Prepare and analyze current window
        WindowData windowData = prevOverlapData;
        analyzeWindow(window, windowStart, windowData, nextOverlapData);

        // Update windowData
        windowData.windowStart = windowStart + absPos;
        windowData.currentWindowSize = currentWindowSize;
        windows.emplace_back(windowData); // Add to the vector of windows

        // Pass and reset overlap data
        prevOverlapData = nextOverlapData;
        nextOverlapData = WindowData(); // Reset nextOverlapData for the next iteration

        // Prepare next window
        windowStart += step;
        if (windowStart >= sequence.size()) {
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

    return windows;
}


void Teloscope::writeBEDFile(std::ofstream& shannonFile, std::ofstream& gcContentFile,
                std::unordered_map<std::string, std::ofstream>& patternMatchFiles,
                std::unordered_map<std::string, std::ofstream>& patternCountFiles,
                std::unordered_map<std::string, std::ofstream>& patternDensityFiles) {
    
    if (!userInput.keepWindowData) { // If windowData is not stored, return
        return;
    }

    for (const auto& windowData : allWindows) {
        unsigned int seqPos;
        std::string header;
        std::vector<WindowData> windows;
        std::tie(seqPos, header, windows) = windowData; // Unpack the tuple

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
                    for (auto pos : data.wMatches) {
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
    std::ofstream telomereBEDFile; // CHECK
    std::ofstream telomereCountFile; // CHECK
    std::unordered_map<std::string, std::ofstream> patternMatchFiles; // Jack: replace with vector to reduce cache locality?
    std::unordered_map<std::string, std::ofstream> patternCountFiles;
    std::unordered_map<std::string, std::ofstream> patternDensityFiles;
    std::cout << "Reporting window matches and metrics in BED/BEDgraphs...\n";

    // Open files once if their modes are enabled
    if (userInput.modeEntropy) {
        shannonFile.open(outRoute + "/shannonEntropy.bedgraph");
    }

    if (userInput.modeGC) {
        gcContentFile.open(outRoute + "/gcContent.bedgraph");
    }

    if (userInput.modeMatch) {
        // telomereBEDFile.open(outRoute + "/telomere_blocks.bed"); // CHECK
        // telomereCountFile.open(outRoute + "/telomere_block_counts.txt"); // CHECK

        for (const auto& pattern : userInput.patterns) {
            patternMatchFiles[pattern].open(outRoute + "/" + pattern + "_matches.bed");
            patternCountFiles[pattern].open(outRoute + "/" + pattern + "_count.bedgraph");
            patternDensityFiles[pattern].open(outRoute + "/" + pattern + "_density.bedgraph");
        }
    }

    // Write data for each window
    writeBEDFile(shannonFile, gcContentFile, 
                // telomereBEDFile, telomereCountFile, // CHECK
                patternMatchFiles, patternCountFiles, patternDensityFiles);

    // Close all files once
    if (userInput.modeEntropy) {
        shannonFile.close();
    }
    if (userInput.modeGC) {
        gcContentFile.close();
    }
    if (userInput.modeMatch) {
        // telomereBEDFile.close(); // CHECK
        // telomereCountFile.close(); // CHECK

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