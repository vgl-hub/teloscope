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

    for (uint64_t i = 0; i < window.size(); ++i) {
        windowData.nucleotideCounts[window[i]]++; // For GC/entropy

        auto current = trie.getRoot();
        
        for (uint64_t j = i; j < window.size(); ++j) {
            
            if (!trie.hasChild(current, window[j])) break;
            current = trie.getChild(current, window[j]);

            if (current->isEndOfWord) {
                std::string pattern = window.substr(i, j - i + 1);
                windowData.patternCounts[pattern]++; // Count all matches

                if (userInput.windowSize == userInput.step || windowStart == 0 || j >= userInput.windowSize - userInput.step) {
                    windowData.patternBEDData.emplace_back(windowStart + i, pattern); // Correct by absPos here?
                }

                // std::vector<uint64_t>>& patternIndexMatches
                // patternIndexMatches[pattern].push_back(i);
                // std::unordered_map<std::string, std::vector<uint64_t>> &patternIndexMatches
                // patternIndexMatches[pattern].push_back(windowStart + i);
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

    while (windowStart < sequence.size()) {
        currentWindowSize = std::min(userInput.windowSize, static_cast<uint32_t>(sequence.size() - windowStart));
        
        WindowData windowData;
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




// void analyzeSegment(std::string header, std::string &sequence, UserInputTeloscope userInput, uint64_t absPos, unsigned int pathId) {

//     std::vector<std::tuple<uint64_t, float>> GCData, entropyData;
//     std::vector<std::tuple<uint64_t, std::string>> patternBEDData;
//     std::map<std::string, std::vector<std::tuple<uint64_t, uint64_t>>> patternCountData;
//     std::map<std::string, std::vector<std::tuple<uint64_t, float>>> patternDensityData;

//     std::unordered_map<char, uint64_t> nucleotideCounts;
//     std::unordered_map<std::string, uint32_t> patternCounts;

//     std::string window = sequence.substr(0, userInput.windowSize);
//     uint64_t windowStart = 0;
//     uint32_t currentWindowSize;

//     while (windowStart < sequence.size()) { // While we are not in the last window
//         currentWindowSize = (userInput.windowSize <= sequence.size() - windowStart) ? userInput.windowSize : (sequence.size() - windowStart);

//         patternCounts.clear();
//         analyzeWindow(window, windowStart, patternBEDData, userInput, nucleotideCounts, patternCounts); //this has to return a populated container

//         for (const auto& [pattern, count] : patternCounts) {
//             float density = static_cast<float>(count * pattern.size()) / currentWindowSize;
//             patternDensityData[pattern].emplace_back(windowStart + absPos, density);
//             patternCountData[pattern].emplace_back(windowStart + absPos, count);
//         }

//         GCData.emplace_back(windowStart + absPos, getGCContent(nucleotideCounts, currentWindowSize));
//         entropyData.emplace_back(windowStart + absPos, getShannonEntropy(nucleotideCounts, currentWindowSize));

//         windowStart += userInput.step;
//         if (currentWindowSize == userInput.windowSize) {
//             window = window.substr(userInput.step) + sequence.substr(windowStart + userInput.windowSize - userInput.step, userInput.step);
//         } else {
//             break;
//         }
//     }
    

// }




// void analyzeSegment(const std::string& sequence, uint64_t absPos, const UserInputTeloscope& userInput,
//                     std::shared_ptr<TrieNode> trieRoot, PathData& pathData) { // PathData is a struc

//     // PathData pathData(pathId); // recycle pathProperties
//     //  PathData pathData; // Initialize local pathData structure. Decide which is better
//     std::map<std::string, std::vector<uint64_t>> patternIndexMatches;
//     std::map<std::string, uint64_t> patternCounts;
    
//     // Copilot #1: Sliding windows
//     for (size_t start = 0; start < sequence.size(); start += userInput.step) {
//         size_t end = std::min(start + userInput.windowSize, sequence.size());
//         std::string window = sequence.substr(start, end - start);
        
//         analyzeWindow(window, trieRoot, absPos + start, patternIndexMatches, patternCounts); // do we need to pass pathData
        
//         // Calculate metrics for this window
//         auto nucleotideCounts = calculateNucleotideCounts(window);  // we don't need this as is calculated already
//         float gcContent = calculateGCContent(nucleotideCounts, window.size());
//         float shannonEntropy = calculateShannonEntropy(nucleotideCounts, window.size());
        
//         // Update PathData with window analysis
//         pathData.updateWindowData(start + absPos, gcContent, shannonEntropy, patternIndexMatches, patternCounts);
//     }

//     // Copilot #2: Sliding windows. Using void analyzeSegment(std::string header, std::string& sequence, UserInputTeloscope userInput, uint64_t absPos, unsigned int pathId) {
//     for (size_t windowStart = 0; windowStart < sequence.length(); windowStart += userInput.step) {
//         WindowData windowData; // Conceptual structure for window-specific data
//         std::string window = sequence.substr(windowStart, userInput.windowSize);
        
//         // Populate windowData
//         analyzeWindow(window, windowStart + absPos, root, windowData); // This function embodies GC content, entropy calculation, and pattern matching
        
//         // Add windowData to pathData
//         pathData.addWindowData(windowData);
//     }

//     // Add pathData to a global or thread-safe container
//     // Note: Thread safety considerations deferred
//     assemblyData.addPathData(pathData);
// }

// // Generate BED files based on consolidated pathData
// generateBEDFiles(assemblyData);


//// Last stable w/o OOP
// void analyzeWindow(std::shared_ptr<TrieNode> root, const std::string &window, uint64_t windowStart, 
//                     std::vector<std::tuple<uint64_t, std::string>> &patternBEDData, const UserInputTeloscope& userInput,
//                     std::unordered_map<char, uint64_t> &nucleotideCounts, std::unordered_map<std::string, uint32_t> &patternCounts) {

//     nucleotideCounts = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};

//     for (uint64_t i = 0; i < window.size(); ++i) {
//         nucleotideCounts[window[i]]++; // For GC/entropy

//         auto current = root;
//         for (uint64_t j = i; j < window.size(); ++j) {

//             if (current->children.find(window[j]) == current->children.end()) break;
//             current = current->children[window[j]];

//             if (current->isEndOfWord) {
//                 std::string pattern = window.substr(i, j - i + 1);
//                 patternCounts[pattern]++; // Count all matches

//                 if (userInput.windowSize == userInput.step || windowStart == 0 || j >= userInput.windowSize - userInput.step) {
//                     patternBEDData.emplace_back(windowStart + i, pattern);
//                 }

//                 // std::vector<uint64_t>>& patternIndexMatches
//                 // patternIndexMatches[pattern].push_back(i);
//                 // std::unordered_map<std::string, std::vector<uint64_t>> &patternIndexMatches
//                 // patternIndexMatches[pattern].push_back(windowStart + i);

//             }
//         }
//     }
// }

//// Last stable w/o OOP
// template <typename T>
// void BEDFileGenerator::generateBEDFile(const std::vector<std::tuple<uint64_t, T>>& data, const std::string& fileName, const std::string& sequence) {
//     std::string cleanedHeader = cleanString(this->header); // Jack: Put directly
//     std::string bedFileName = outRoute + "/teloscope_" + fileName + (typeid(T) == typeid(std::string) ? ".bed" : ".bedgraph");
//     std::ofstream bedFile(bedFileName, std::ios::out | std::ios::app); // Open file in append mode

//     if (!bedFile.is_open()) {
//         std::cerr << "Failed to open file: " << bedFileName << '\n';
//         return;
//     }

//     for (const auto& entry : data) {
//         uint64_t start = std::get<0>(entry);
//         T value = std::get<1>(entry);
//         uint64_t end;

//         if constexpr (std::is_same<T, std::string>::value) {
//             end = start + value.size() - 1; // For patterns, use the string size
//         } else {
//             end = (start + this->userInput.windowSize - 1 < sequence.size()) ? start + this->userInput.windowSize - 1 : sequence.size() - 1;
//         }

//         bedFile << cleanedHeader << "\t" << start + this->absPos << "\t" << end + this->absPos << "\t" << value << "\n";

//     }

//     bedFile.close();
// }

// //// Last stable w/o OOP 
// void analyzeSegment(std::string header, std::string &sequence, UserInputTeloscope userInput, uint64_t absPos) {

//     Trie trie;
//     for (const auto& pattern : userInput.patterns) {
//         trie.insertPattern(pattern);
//     }
//     std::vector<std::tuple<uint64_t, float>> GCData, entropyData;
//     std::vector<std::tuple<uint64_t, std::string>> patternBEDData;
//     std::map<std::string, std::vector<std::tuple<uint64_t, uint64_t>>> patternCountData;
//     std::map<std::string, std::vector<std::tuple<uint64_t, float>>> patternDensityData;

//     std::map<char, uint64_t> nucleotideCounts;
//     std::unordered_map<std::string, uint32_t> patternCounts;

//     std::string window = sequence.substr(0, userInput.windowSize);
//     uint64_t windowStart = 0;
//     uint32_t currentWindowSize;

//     while (windowStart < sequence.size()) {
//         currentWindowSize = (userInput.windowSize <= sequence.size() - windowStart) ? userInput.windowSize : (sequence.size() - windowStart);

//         patternCounts.clear();
//         analyzeWindow(window, windowStart, patternBEDData, userInput, nucleotideCounts, patternCounts);

//         for (const auto& [pattern, count] : patternCounts) {
//             float density = static_cast<float>(count * pattern.size()) / currentWindowSize;
//             patternDensityData[pattern].emplace_back(windowStart + absPos, density);
//             patternCountData[pattern].emplace_back(windowStart + absPos, count);
//         }

//         GCData.emplace_back(windowStart + absPos, getGCContent(nucleotideCounts, currentWindowSize));
//         entropyData.emplace_back(windowStart + absPos, getShannonEntropy(nucleotideCounts, currentWindowSize));

//         windowStart += userInput.step;
//         if (currentWindowSize == userInput.windowSize) {
//             window = window.substr(userInput.step) + sequence.substr(windowStart + userInput.windowSize - userInput.step, userInput.step);
//         } else {
//             break;
//         }
//     }

//     BEDFileGenerator bedFileGenerator(header, userInput, absPos);  // Instantiate BEDFileGenerator

//     bedFileGenerator.generateBEDFile(GCData, "window_gc", sequence);
//     bedFileGenerator.generateBEDFile(entropyData, "window_entropy", sequence);
//     bedFileGenerator.generateBEDFile(patternBEDData, "pattern_matching", sequence);

//     for (const auto& [pattern, countData] : patternCountData) {
//         bedFileGenerator.generateBEDFile(countData, pattern + "_count", sequence);
//     }
//     for (const auto& [pattern, densityData] : patternDensityData) {
//         bedFileGenerator.generateBEDFile(densityData, pattern + "_density", sequence);
//     }
// }
