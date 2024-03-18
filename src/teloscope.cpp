#include <iostream>
#include <fstream>
#include <sstream> // check
#include <map>
#include <stdint.h> // what's this for?
#include <vector>
#include <algorithm> // cleanString
#include <array> // check
#include <cmath>
#include <type_traits> // generateBEDFile
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
}


void Trie::findPatternsInWindow(const std::string &window, uint64_t windowStart,
                                std::vector<std::tuple<uint64_t, std::string>> &patternBEDData, const UserInputTeloscope& userInput,
                                std::map<char, uint64_t> &nucleotideCounts, std::unordered_map<std::string, uint32_t> &patternCounts) {

    nucleotideCounts = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};

    for (uint64_t i = 0; i < window.size(); ++i) {
        nucleotideCounts[window[i]]++; // For GC/entropy

        auto current = root;
        for (uint64_t j = i; j < window.size(); ++j) {

            if (current->children.find(window[j]) == current->children.end()) break;
            current = current->children[window[j]];

            if (current->isEndOfWord) {
                std::string pattern = window.substr(i, j - i + 1);
                patternCounts[pattern]++; // Count all matches

                if (userInput.windowSize == userInput.step || windowStart == 0 || j >= userInput.windowSize - userInput.step) {
                    patternBEDData.emplace_back(windowStart + i, pattern);
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
                    const std::string& fileName, const UserInputTeloscope& userInput, std::string &sequence) {

    std::string cleanedHeader = cleanString(header);
    std::string bedFileName = outRoute + "/teloscope_" + fileName + (typeid(T) == typeid(std::string) ? ".bed" : ".bedgraph");
    std::ofstream bedFile(bedFileName, std::ios::out | std::ios::app); // Open file in append mode
    
    if (!bedFile.is_open()) {
        std::cerr << "Failed to open file: " << bedFileName << '\n';
        return;
    }

    for (const auto& entry : data) {
        uint64_t start = std::get<0>(entry);
        T value = std::get<1>(entry);
        uint64_t end;

        if constexpr (std::is_same<T, std::string>::value) {
            end = start + value.size() - 1; // For patterns, use the string size
        } else {
            end = (start + userInput.windowSize - 1 < sequence.size()) ? start + userInput.windowSize - 1 : sequence.size() - 1; // Otherwise, use the window size
        }

        bedFile << cleanedHeader << "\t" << start + userInput.absPos << "\t" << end + userInput.absPos<< "\t" << value << "\n";
    }

    bedFile.close();
}


void findTelomeres(std::string header, std::string &sequence, UserInputTeloscope userInput) { // const UserInputTeloscope& userInput
   
    Trie trie;
    for (const auto& pattern : userInput.patterns) {
        trie.insertPattern(pattern);
    }
    std::vector<std::tuple<uint64_t, float>> GCData, entropyData;
    std::vector<std::tuple<uint64_t, std::string>> patternBEDData;
    std::map<std::string, std::vector<std::tuple<uint64_t, uint64_t>>> patternCountData;
    std::map<std::string, std::vector<std::tuple<uint64_t, float>>> patternDensityData;

    std::map<char, uint64_t> nucleotideCounts;
    std::unordered_map<std::string, uint32_t> patternCounts;

    std::string window = sequence.substr(0, userInput.windowSize);
    uint64_t windowStart = 0;
    uint32_t currentWindowSize;

    while (windowStart < sequence.size()) {
        currentWindowSize = (userInput.windowSize <= sequence.size() - windowStart) ? userInput.windowSize : (sequence.size() - windowStart);

        patternCounts.clear();
        trie.findPatternsInWindow(window, windowStart, patternBEDData, userInput, nucleotideCounts, patternCounts);

        for (const auto& [pattern, count] : patternCounts) {
            float density = static_cast<float>(count * pattern.size()) / currentWindowSize;
            patternDensityData[pattern].emplace_back(windowStart, density);
            patternCountData[pattern].emplace_back(windowStart, count);
        }

        GCData.emplace_back(windowStart, getGCContent(nucleotideCounts, currentWindowSize));
        entropyData.emplace_back(windowStart, getShannonEntropy(nucleotideCounts, currentWindowSize));

        windowStart += userInput.step;
        if (currentWindowSize == userInput.windowSize) {
            window = window.substr(userInput.step) + sequence.substr(windowStart + userInput.windowSize - userInput.step, userInput.step);
        } else {
            break;
        }
    }

    generateBEDFile(header, GCData, "window_gc", userInput, sequence);
    generateBEDFile(header, entropyData, "window_entropy", userInput, sequence);
    generateBEDFile(header, patternBEDData, "pattern_matching", userInput, sequence);
    for (const auto& [pattern, countData] : patternCountData) {
        generateBEDFile(header, countData, pattern + "_count", userInput, sequence);
    }
    for (const auto& [pattern, densityData] : patternDensityData) {
        generateBEDFile(header, densityData, pattern + "_density", userInput, sequence);
    }
}