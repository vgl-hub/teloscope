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

#include <parallel_hashmap/phmap.h>

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

BedCoordinates bedCoords;

std::vector<uint64_t> getPatternFrequency(const std::vector<bool>& patternMatches, uint32_t windowSize, uint32_t step) {
    std::vector<uint64_t> patternFreq;
    uint64_t vecBoolLength = patternMatches.size();
    patternFreq.reserve(vecBoolLength);

    if (vecBoolLength < windowSize) {
        return patternFreq;
    }

    uint32_t windowCounts = 0;
    for (uint32_t i = 0; i < windowSize; ++i) {
        if (patternMatches[i]) {
            windowCounts++;
        }
    }
    patternFreq.push_back(windowCounts);

    for (uint64_t i = windowSize; i < vecBoolLength; i += step) {
        windowCounts += patternMatches[i] - patternMatches[i - windowSize];
        patternFreq.push_back(windowCounts);
    }

    return patternFreq;
}


double getShannonEntropy(const std::string& window) { // float, less memory?
    std::array<int, 128> freq = {0};
    for (char nucleotide : window) {
        freq[static_cast<unsigned char>(nucleotide)]++; // giulio: no need to cast, check kmer.h 71
    }
    double entropy = 0.0;
    for (int f : freq) {
        if (f > 0) {
            double probability = static_cast<double>(f) / window.size();
            entropy -= probability * std::log2(probability);
        }
    }
    return entropy;
}

template <typename T>
void generateBEDFile(const std::string& header, const std::vector<std::tuple<uint64_t, uint64_t, T>>& data, const std::string& fileName) {
    std::string bedFileName = "../../output/" + header + "_" + fileName + (std::is_floating_point<T>::value ? ".bedgraph" : ".bed");
    std::ofstream bedFile(bedFileName, std::ios::out);
    if (!bedFile.is_open()) {
        std::cerr << "Failed to open file: " << bedFileName << std::endl;
        return;
    }

    for (const auto& entry : data) {
        uint64_t start = std::get<0>(entry);
        uint64_t end = std::get<1>(entry);
        T value = std::get<2>(entry);

        bedFile << header << "\t" << start << "\t" << end << "\t" << value << std::endl;
    }

    bedFile.close();
}

void findTelomeres(std::string header, std::string &sequence, UserInputTeloscope userInput) {
    uint64_t segLength = sequence.size();
    uint32_t windowSize = static_cast<uint32_t>(userInput.windowSize);
    uint32_t step = static_cast<uint32_t>(userInput.step);

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> patternBEDData;
    std::vector<std::tuple<uint64_t, uint64_t, double>> entropyData;
    std::map<std::string, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>> patternCountData;
    std::map<std::string, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>> patternFractionData;

    for (uint64_t windowStart = 0; windowStart <= segLength - windowSize; windowStart += step) {
        std::string window = sequence.substr(windowStart, windowSize);
        double entropy = getShannonEntropy(window);
        std::cout << "\nShannon Entropy for window [" << windowStart << ", " << (windowStart + windowSize - 1) << "]: " << entropy << std::endl;
        entropyData.emplace_back(windowStart, windowStart + windowSize - 1, entropy);
        // double totalPatternFraction = 0.0;

        std::map<std::string, uint64_t> windowPatternCounts;
        for (const auto& pattern : userInput.patterns) {
            uint64_t patternLength = pattern.size();
            uint64_t patternCount = 0;

            // Diagnostic Print
            std::cout << "Analyzing Pattern: " << pattern << std::endl;
            
            // Merged loop
            for (uint64_t i = windowStart; i < windowStart + windowSize && i + patternLength <= segLength; ++i) {
                if (pattern == sequence.substr(i, patternLength) && i + patternLength <= windowStart + windowSize) {

                    // Diagnostic Print
                    std::cout << "Pattern Match at Position: " << i << std::endl;

                    patternCount++;
                    patternBEDData.emplace_back(i, i + patternLength - 1, pattern);
                    i += patternLength - 1; // Skip ahead to avoid duplicate entries
                }
            }
            windowPatternCounts[pattern] = patternCount;

            // Print pattern fraction for the window
            double patternFraction = static_cast<double>(patternCount * patternLength) / windowSize;
            std::cout << "Window " << windowStart << ": Pattern \"" << pattern << "\" = " << patternFraction << std::endl;
        }

        // Arrange pattern count data for each window
        for (const auto& [pattern, count] : windowPatternCounts) {
            patternCountData[pattern].emplace_back(windowStart, windowStart + windowSize - 1, count);
        }
    }

    // Generate BED and BEDgraph files
    generateBEDFile(header, patternBEDData, "pattern_mapping");
    generateBEDFile(header, entropyData, "window_entropy");
    for (const auto& [pattern, data] : patternCountData) {
        generateBEDFile(header, data, pattern + "_count");
    }
}

// test
int main() {
    std::string header = "H1.scaffold_401";
    std::string sequence = "TTAGGGGCCCCTTCGATCGACCCTTTAGGG";
    UserInputTeloscope userInput;
    userInput.patterns = {"TTAGGG", "CCCT"};
    userInput.windowSize = 10;
    userInput.step = 5;

    findTelomeres(header, sequence, userInput);

    return 0;
}