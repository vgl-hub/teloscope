#include <iostream>
#include <fstream>
#include <sstream> // check
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

    if (pattern.size() > longestPatternSize) {
        longestPatternSize = pattern.size();
    }
}


std::string cleanString(const std::string& input) {
    std::string output = input;
    output.erase(std::remove(output.begin(), output.end(), '\r'), output.end());
    return output;
}


float Teloscope::getShannonEntropy(const std::unordered_map<char, uint64_t>& nucleotideCounts, uint32_t windowSize) {
    float entropy = 0.0;
    for (auto &[nucleotide, count] : nucleotideCounts) {
        if (count > 0) {
            float probability = static_cast<float>(count) / windowSize;
            entropy -= probability * std::log2(probability);
        }
    }
    return entropy;
}


float Teloscope::getGCContent(const std::unordered_map<char, uint64_t>& nucleotideCounts, uint32_t windowSize) {
    uint64_t gcCount = nucleotideCounts.at('G') + nucleotideCounts.at('C');
    return float(gcCount) / windowSize * 100.0;
}


void Teloscope::analyzeWindow(const std::string &window, uint64_t windowStart, WindowData& windowData) {

    windowData.windowStart = windowStart;
    unsigned short int longestPatternSize = this->trie.getLongestPatternSize(); // to implement

    for (uint64_t i = 0; i < window.size(); ++i) { // For each nucleotide in the window
        windowData.nucleotideCounts[window[i]]++; // For GC/entropy

        auto current = trie.getRoot();
        uint64_t scanLimit = std::min(i + longestPatternSize, window.size());
        
        for (uint64_t j = i; j < scanLimit; ++j) { // Only scan positions in range of patterns
            
            if (!trie.hasChild(current, window[j])) break;  
            current = trie.getChild(current, window[j]);

            if (current->isEndOfWord) {
                std::string pattern = window.substr(i, j - i + 1);
                windowData.patternCounts[pattern]++; // Count all matches

                if (userInput.windowSize == userInput.step || windowStart == 0 || j >= userInput.windowSize - userInput.step) {
                    windowData.patternBEDData.emplace_back(windowStart + i, pattern); // Correct by absPos here?
                }


            }
        }
    }

    windowData.gcContent = getGCContent(windowData.nucleotideCounts, window.size());
    windowData.shannonEntropy = getShannonEntropy(windowData.nucleotideCounts, window.size());
}


// std::vector<WindowData> Teloscope::analyzeSegment(std::string header, std::string &sequence, UserInputTeloscope userInput, uint64_t absPos, unsigned int pathId) {
std::vector<WindowData> Teloscope::analyzeSegment(std::string &sequence, UserInputTeloscope userInput, uint64_t absPos) {

    uint64_t windowStart = 0;
    uint32_t currentWindowSize;
    std::string window = sequence.substr(0, userInput.windowSize);
    std::vector<WindowData> windows;

    // for (uint64_t windowStart = 0; windowStart + userInput.windowSize <= sequence.size(); windowStart += userInput.step) {
    //     WindowData windowData;
    //     windowData.windowStart = windowStart + absPos;

    //     analyzeWindow(window, windowStart, windowData);
        

    while (windowStart < sequence.size()) {
        currentWindowSize = std::min(userInput.windowSize, static_cast<uint32_t>(sequence.size() - windowStart));
        
        WindowData windowData; // create instance per window
        windowData.patternCounts.clear(); // otherwise it accumulates? 
        // analyzeWindow(window, windowStart + absPos, windowData);
        analyzeWindow(window, windowStart, windowData);
        windowData.windowStart = windowStart + absPos;

        for (const auto& [pattern, count] : windowData.patternCounts) {
            float density = static_cast<float>(count * pattern.size()) / currentWindowSize;
            windowData.patternDensityData[pattern].emplace_back(windowStart + absPos, density);
            windowData.patternCountData[pattern].emplace_back(windowStart + absPos, count);
        }

        windows.push_back(windowData);

        windowStart += userInput.step;
        if (currentWindowSize == userInput.windowSize) {
            window = window.substr(userInput.step) + sequence.substr(windowStart + userInput.windowSize - userInput.step, userInput.step);
        } else {
            break;
        }
    }
    
    return windows;
    
}