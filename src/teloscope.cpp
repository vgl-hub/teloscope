#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <array>
#include <cmath>

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
        freq[static_cast<unsigned char>(nucleotide)]++;
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

void findTelomeres(std::string header, std::string &sequence, UserInputTeloscope userInput) {
    uint64_t segLength = sequence.size();
    uint32_t windowSize = static_cast<uint32_t>(userInput.windowSize);
    uint32_t step = static_cast<uint32_t>(userInput.step);

    // Map of pattern matches
    std::map<std::string, std::vector<bool>> patternMatchesMap;
    uint64_t vecBoolLength = 0; 
    for (const auto& pattern : userInput.patterns) {
        vecBoolLength = segLength - pattern.size() + 1; 
        patternMatchesMap[pattern] = std::vector<bool>(vecBoolLength, false);
        for (uint64_t i = 0; i <= segLength - pattern.size(); ++i) {
            if (pattern == sequence.substr(i, pattern.size())) {
                patternMatchesMap[pattern][i] = true;
            }
        }
    }

    // Diagnostic Print
    std::cout << "Pattern Matches Map Size: " << vecBoolLength << std::endl;

    // BED File Setup
    std::string bedFileName = "../../output/" + header + "_all_patterns_matches.bed";
    std::ofstream bedFile(bedFileName, std::ios::out); // Open in write mode
    if (!bedFile.is_open()) {
        std::cerr << "Failed to open BED file: " << bedFileName << std::endl;
        return;
    }

    // Looping over the sliding window
    for (uint64_t windowStart = 0; windowStart <= segLength - windowSize; windowStart += step) {
        std::string window = sequence.substr(windowStart, windowSize);
        double entropy = getShannonEntropy(window);
        std::cout << "\nShannon Entropy for window [" << windowStart + 1 << ", " << (windowStart + windowSize) << "]: " << entropy << std::endl;

        double totalPatternFraction = 0.0;
        for (const auto& pattern : userInput.patterns) {
            uint64_t patternLength = pattern.size();
            std::vector<bool> currentWindowMatches(windowSize, false);

            // Diagnostic Print
            std::cout << "Analyzing Pattern: " << pattern << std::endl;            

            for (uint64_t i = windowStart; i < windowStart + windowSize && i < vecBoolLength; ++i) {
                if (patternMatchesMap[pattern][i]) {
                    currentWindowMatches[i - windowStart] = true;
                    
                    // Diagnostic Print
                    std::cout << "Pattern Match at Position: " << i << std::endl;
                }
            }
            uint64_t patternCount = std::count(currentWindowMatches.begin(), currentWindowMatches.end(), true);
            double patternFraction = static_cast<double>(patternCount * patternLength) / windowSize;
            totalPatternFraction += patternFraction;

            std::cout << "Window " << windowStart + 1 << ": Pattern \"" << pattern << "\" = " << patternFraction << std::endl;

            // BED file logic here
            for (uint64_t i = windowStart; i < windowStart + windowSize && i < vecBoolLength; ++i) {
                if (patternMatchesMap[pattern][i]) {
                    bedFile << header << "\t" << i << "\t" << i + patternLength << "\t" << pattern << std::endl;
                    i += patternLength - 1; // Skip ahead to avoid duplicate entries
                }
            }
        }

        double noneFraction = std::max(1.0 - totalPatternFraction, 0.0); // Ensure non-negative
        std::cout << "Window " << windowStart + 1 << ": None = " << noneFraction << std::endl;
    }

    bedFile.close();
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