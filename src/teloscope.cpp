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

// UserInputTeloscope userInput;

float getShannonEntropy(const std::string& window) {
    std::array<int, 128> freq = {0};
    for (char nucleotide : window) {
        freq[(nucleotide)]++; // giulio: no need to cast, check kmer.h 71
    }
    float entropy = 0.0;
    for (int f : freq) {
        if (f > 0) {
            float probability = static_cast<float>(f) / window.size();
            entropy -= probability * std::log2(probability);
        }
    }
    return entropy;
}

float getGCContent(const std::string& window) {
    uint64_t gcCount = 0;

    gcCount = std::count_if(window.begin(), window.end(), [](char n) { return n == 'G' || n == 'C'; });

    return float(gcCount) / window.size() * 100.0; 
}

template <typename T>
void generateBEDFile(const std::string& header, const std::vector<std::tuple<uint64_t, uint64_t, T>>& data, const std::string& fileName) {
    std::string bedFileName = "../../output/" + header + "_" + fileName + (typeid(T) == typeid(std::string) ? ".bed" : ".bedgraph");
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
    uint32_t windowSize = userInput.windowSize;
    uint32_t step = userInput.step;

    std::vector<std::tuple<uint64_t, uint64_t, float>> GCData;
    std::vector<std::tuple<uint64_t, uint64_t, float>> entropyData;
    std::vector<std::tuple<uint64_t, uint64_t, std::string>> patternBEDData;
    std::map<std::string, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>> patternCountData;
    std::map<std::string, std::vector<std::tuple<uint64_t, uint64_t, float>>> patternFractionData;

    std::string window = sequence.substr(0, windowSize);
    uint64_t windowStart = 0;

    // for (uint64_t windowStart = 0; windowStart <= segLength - windowSize; windowStart += step) {
    while (windowStart + windowSize <= segLength) {

        float GC = getGCContent(window);
        GCData.emplace_back(windowStart, windowStart + windowSize - 1, GC);

        float entropy = getShannonEntropy(window);
        entropyData.emplace_back(windowStart, windowStart + windowSize - 1, entropy);

        for (const auto& pattern : userInput.patterns) {
            uint64_t patternLength = pattern.size();
            uint64_t patternCount = 0;

            for (uint64_t i = 0; i + patternLength <= windowSize; ++i) {
                if (pattern == window.substr(i, patternLength)) {
                    patternCount++;
                    patternBEDData.emplace_back(windowStart + i, windowStart + i + patternLength - 1, pattern);
                }
            }

            float patternFraction = static_cast<float>(patternCount * patternLength) / windowSize;
            patternFractionData[pattern].emplace_back(windowStart, windowStart + windowSize - 1, patternFraction);
            patternCountData[pattern].emplace_back(windowStart, windowStart + windowSize - 1, patternCount);
        }

        // Update window for next iteration
        windowStart += step;
        if (windowStart + windowSize <= segLength) {
            window = window.substr(step) + sequence.substr(windowStart + windowSize - step, step);
        }
    }

    // Generate BED and BEDgraph files
    
    generateBEDFile(header, GCData, "window_gc");
    generateBEDFile(header, entropyData, "window_entropy");
    generateBEDFile(header, patternBEDData, "pattern_mapping");
    for (const auto& [pattern, countData] : patternCountData) {
        generateBEDFile(header, countData, pattern + "_count");
    }
    for (const auto& [pattern, fractionData] : patternFractionData) {
        generateBEDFile(header, fractionData, pattern + "_fraction");
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

    auto start = std::chrono::high_resolution_clock::now();
    findTelomeres(header, sequence, userInput);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " milliseconds" << '\n';

    return 0;
}