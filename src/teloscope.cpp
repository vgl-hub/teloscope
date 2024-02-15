#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>
#include <chrono>
#include <memory>
#include <unordered_map>

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

struct TrieNode {
    std::unordered_map<char, std::shared_ptr<TrieNode>> children;
    bool isEndOfWord = false;
};

void insertPattern(std::shared_ptr<TrieNode> root, const std::string &pattern) {
    auto current = root;
    for (char ch : pattern) {
        if (current->children.find(ch) == current->children.end()) {
            current->children[ch] = std::make_shared<TrieNode>();
        }
        current = current->children[ch];
    }
    current->isEndOfWord = true;
}

void findPatternsInWindow(std::shared_ptr<TrieNode> root, const std::string &window,
                        uint64_t windowStart, std::vector<std::tuple<uint64_t, std::string>> &patternBEDData,
                        std::map<std::string, uint64_t> &lastPatternPositions, uint32_t step,
                        std::map<char, uint64_t> &nucleotideCounts, std::unordered_map<std::string, uint32_t> &patternCounts) {
    
    nucleotideCounts = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
    
    for (uint64_t i = 0; i < window.size(); ++i) {
        nucleotideCounts[window[i]]++;

        auto current = root;
        for (uint64_t j = i; j < window.size(); ++j) {

            if (current->children.find(window[j]) == current->children.end()) break;
            current = current->children[window[j]];

            if (current->isEndOfWord) {
                std::string pattern = window.substr(i, j - i + 1);
                patternCounts[pattern]++; // outside to count them per window

                if (lastPatternPositions[pattern] <= windowStart + i) {
                    patternBEDData.emplace_back(windowStart + i, pattern);
                    lastPatternPositions[pattern] = windowStart + i + step - 1;
                }
            }
        }
    }
}



float getShannonEntropy(const std::map<char, uint64_t>& nucleotideCounts, uint32_t windowSize) {
    float entropy = 0.0;
    for (auto &[nucleotide, count] : nucleotideCounts) {
        if (count > 0) {
            float probability = static_cast<float>(count) / windowSize;
            entropy -= probability * std::log2(probability);
        }
    }
    return entropy;
}

float getGCContent(const std::map<char, uint64_t>& nucleotideCounts, uint32_t windowSize) {
    uint64_t gcCount = nucleotideCounts.at('G') + nucleotideCounts.at('C');
    return float(gcCount) / windowSize * 100.0;
}

std::string cleanString(const std::string& input) {
    std::string output = input;
    output.erase(std::remove(output.begin(), output.end(), '\r'), output.end());
    return output;
}


template <typename T>
void generateBEDFile(const std::string& header, const std::vector<std::tuple<uint64_t, T>>& data, 
                    const std::string& fileName, uint32_t windowSize, uint64_t segLength) {    
    std::string cleanedHeader = cleanString(header); // Clean the header string
    std::string bedFileName = "../../output/" + cleanedHeader + "_" + fileName + (typeid(T) == typeid(std::string) ? ".bed" : ".bedgraph");
    std::ofstream bedFile(bedFileName, std::ios::out);
    
    if (!bedFile.is_open()) {
        std::cerr << "Failed to open file: " << bedFileName << '\n';
        return;
    }

    for (const auto& entry : data) {
        uint64_t start = std::get<0>(entry);
        T value = std::get<1>(entry);
        uint64_t end;

        if constexpr (std::is_same<T, std::string>::value) {
            end = start + value.size() - 1; // For pattern data, use the string size
        } else {
            end = (start + windowSize - 1 < segLength) ? start + windowSize - 1 : segLength - 1;; // For other data, use the window size
        }

        bedFile << cleanedHeader << "\t" << start << "\t" << end << "\t" << value << "\n";
    }

    bedFile.close();
}

void findTelomeres(std::string header, std::string &sequence, UserInputTeloscope userInput) {
    uint64_t segLength = sequence.size();
    uint32_t windowSize = userInput.windowSize;
    uint32_t step = userInput.step;

    auto root = std::make_shared<TrieNode>();
    for (const auto& pattern : userInput.patterns) {
        insertPattern(root, pattern);
    }
    
    std::vector<std::tuple<uint64_t, float>> GCData, entropyData;
    std::vector<std::tuple<uint64_t, std::string>> patternBEDData;
    std::map<std::string, std::vector<std::tuple<uint64_t, uint64_t>>> patternCountData;

    std::map<std::string, uint64_t> lastPatternPositions;
    std::map<char, uint64_t> nucleotideCounts;
    std::unordered_map<std::string, uint32_t> patternCounts;

    std::string window = sequence.substr(0, windowSize);
    uint64_t windowStart = 0;
    uint32_t currentWindowSize = userInput.windowSize;

    while (windowStart < segLength) {
        currentWindowSize = (windowSize <= segLength - windowStart) ? windowSize : (segLength - windowStart);

        patternCounts.clear();
        findPatternsInWindow(root, window, windowStart, patternBEDData, lastPatternPositions, step, nucleotideCounts, patternCounts);

        for (const auto& [pattern, count] : patternCounts) {
            patternCountData[pattern].emplace_back(windowStart, count);
        }

        GCData.emplace_back(windowStart, getGCContent(nucleotideCounts, currentWindowSize));
        entropyData.emplace_back(windowStart, getShannonEntropy(nucleotideCounts, currentWindowSize));
        
        windowStart += step;
        if (currentWindowSize == windowSize) {
            window = window.substr(step) + sequence.substr(windowStart + windowSize - step, step);
        } else {
            break;
        }
    }

    // Example call within findTelomeres
    generateBEDFile(header, GCData, "window_gc", windowSize, segLength);
    generateBEDFile(header, entropyData, "window_entropy", windowSize, segLength);
    generateBEDFile(header, patternBEDData, "pattern_mapping", windowSize, segLength);
    for (const auto& [pattern, countData] : patternCountData) {
        generateBEDFile(header, countData, pattern + "_count", windowSize, segLength);
    }
}