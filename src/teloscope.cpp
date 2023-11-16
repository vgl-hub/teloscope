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

// #include "bed.h"
#include "gfalibs/include/bed.h"
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
    uint64_t seqLength = patternMatches.size();
    patternFreq.reserve(seqLength);

    if (seqLength < windowSize) {
        return patternFreq;
    }

    uint32_t windowCounts = 0;
    for (uint32_t i = 0; i < windowSize; ++i) {
        if (patternMatches[i]) {
            windowCounts++;
        }
    }
    patternFreq.push_back(windowCounts);

    for (uint64_t i = windowSize; i < seqLength; i += step) {
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

    //test
    std::string dummySequence = "TTAGGGCCCCTTTAGGGTTAGGGCCCCTTTAGGGTTAGGGCCCCTTTAGGG"; 
    header = "test"; 
    sequence = dummySequence; 
    userInput.patterns = {"TTAGGG", "CCCCT"}; 
    userInput.windowSize = 10; 
    userInput.step = 5;


    uint64_t seqLength = sequence.size();
    uint32_t windowSize = static_cast<uint32_t>(userInput.windowSize);
    uint32_t step = static_cast<uint32_t>(userInput.step);
    

    // Map of pattern matches
    std::map<std::string, std::vector<bool>> patternMatchesMap;
    for (const auto& pattern : userInput.patterns) {
        patternMatchesMap[pattern] = std::vector<bool>(seqLength - pattern.size() + 1, false);
    }

    // Single pass over the sequence
    for (uint64_t i = 0; i <= seqLength; ++i) {
        for (const auto& pattern : userInput.patterns) {
            if (i + pattern.size() <= seqLength && pattern == sequence.substr(i, pattern.size())) {
                patternMatchesMap[pattern][i] = true; // giulio exit loop when true!!!
            }
        }
    }

    // Pattern range covered in sliding window
    std::vector<uint64_t> totalPatternRange(seqLength, 0);
    for (const auto& pattern : userInput.patterns) {
        std::vector<uint64_t> patternFreq = getPatternFrequency(patternMatchesMap[pattern], windowSize, step);
        for (uint64_t i = 0; i < patternFreq.size(); ++i) {
            totalPatternRange[i] += patternFreq[i] * pattern.size();
        }
    }

    // Calculate and print pattern fractions
    for (const auto& pattern : userInput.patterns) {
        std::vector<uint64_t> patternFreq = getPatternFrequency(patternMatchesMap[pattern], windowSize, step);
        std::cout << "Pattern Frequency Fraction for \"" << pattern << "\": ";

        for (uint64_t i = 0; i < patternFreq.size(); ++i) {
            double patternFraction = static_cast<double>(patternFreq[i] * pattern.size()) / windowSize;
            double noneFraction = 1.0 - static_cast<double>(totalPatternRange[i]) / windowSize;
            std::cout << "Window " << i << ": Pattern = " << patternFraction << ", None = " << noneFraction << " | ";
        }
        std::cout << std::endl;
    }

    for (const auto& pattern : userInput.patterns) {
        std::vector<uint64_t> patternFreq = getPatternFrequency(patternMatchesMap[pattern], windowSize, step);

        // Shannon Entropy
        for (uint64_t i = 0; i < seqLength - windowSize + 1; i += step) {
            std::string window = sequence.substr(i, windowSize);
            double entropy = getShannonEntropy(window);
            std::cout << "Shannon Entropy for window [" << i << ", " << (i + windowSize) << "]: " << entropy << std::endl;
        }

        // Get concatenated BED file for pattern matches
        std::string bedFileName = "./output/" + header + "_all_patterns_matches.bed";
        std::ofstream bedFile(bedFileName, std::ios::app); // Append mode

        for (uint64_t i = 0; i < patternMatchesMap[pattern].size(); ++i) {
            if (patternMatchesMap[pattern][i]) {
                bedFile << header << "\t" << i << "\t" << i + pattern.size() << "\t" << pattern << std::endl;
            }
        }
        bedFile.close();
    }
}