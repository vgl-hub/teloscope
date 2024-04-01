#ifndef TELOSCOPE_H
#define TELOSCOPE_H

#include "input.h" // not in Mac's code
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

std::string cleanString(const std::string& input);

class Trie {
    struct TrieNode {
        std::unordered_map<char, std::shared_ptr<TrieNode>> children;
        bool isEndOfWord = false;
    };

    std::shared_ptr<TrieNode> root;

    unsigned short int longestPatternSize = 0;

public:
    Trie() : root(std::make_shared<TrieNode>()) {}

    void insertPattern(const std::string& pattern);

    std::shared_ptr<TrieNode> getRoot() const {
        return root;
    }

    bool hasChild(const std::shared_ptr<TrieNode>& node, char ch) const { // Giulio: merge with the following method
        return node->children.find(ch) != node->children.end();
    }

    std::shared_ptr<TrieNode> getChild(const std::shared_ptr<TrieNode>& node, char ch) const {
        if (hasChild(node, ch)) {
            return node->children[ch];
        }
        return nullptr;
    }

    // std::shared_ptr<TrieNode> getChildIfExists(const std::shared_ptr<TrieNode>& node, char ch) const {
    //     auto it = node->children.find(ch); // Jack: Test behavior with breaks.
    //     if (it != node->children.end()) {
    //         return it->second;
    //     }
    //     return nullptr;
    // }

    // std::shared_ptr<Trie::TrieNode> Trie::getChildIfExists(const std::shared_ptr<TrieNode>& node, char ch) const {
    //     auto it = node->children.find(ch);
    //     if (it != node->children.end()) {
    //         return it->second; // Directly returns the shared_ptr to the TrieNode if found
    //     }
    //     return nullptr; // Returns nullptr if the character is not a child of the given node
    // }

    unsigned short int getLongestPatternSize() const {
        return longestPatternSize;
    }
};


struct PatternData {
    std::vector<uint64_t> positions; // Match indexes in window
    uint32_t count = 0; // Total pattern count
    float density = 0.0f; // Density of the pattern
};

struct WindowData {
    uint64_t windowStart;
    float gcContent;
    float shannonEntropy;
    std::unordered_map<char, uint64_t> nucleotideCounts;
    std::unordered_map<std::string, PatternData> patternMap; // Condensed pattern data
    
    WindowData() : windowStart(0), gcContent(0.0f), shannonEntropy(0.0f), nucleotideCounts{{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}} {}
};


// struct WindowData {
//     uint64_t windowStart;
//     float gcContent;
//     float shannonEntropy;
//     std::unordered_map<char, uint64_t> nucleotideCounts;
//     std::unordered_map<std::string, uint32_t> patternCounts;
//     std::vector<std::tuple<uint64_t, std::string>> patternBEDData;
//     std::map<std::string, std::vector<std::tuple<uint64_t, uint64_t>>> patternCountData;
//     std::map<std::string, std::vector<std::tuple<uint64_t, float>>> patternDensityData;

    
//     WindowData() : windowStart(0), gcContent(0.0f), shannonEntropy(0.0f) {
//         nucleotideCounts = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
//     }
// };

class Teloscope {

    Trie trie; // Declare trie instance
    UserInputTeloscope userInput; // Declare user input instance
    std::vector<std::pair<unsigned int, std::vector<WindowData>>> allWindows; // Assembly windows

    float getShannonEntropy(const std::unordered_map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);
    float getGCContent(const std::unordered_map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);


public:

    Teloscope(UserInputTeloscope userInput) : userInput(userInput) {
        for (const auto& pattern : userInput.patterns) {
        trie.insertPattern(pattern);
        }
    }

    bool walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);

    void analyzeWindow(const std::string &window, uint64_t windowStart, WindowData& windowData);
    
    // std::vector<WindowData> analyzeSegment(std::string header, std::string &sequence, UserInputTeloscope userInput, uint64_t absPos, unsigned int pathId);
    std::vector<WindowData> analyzeSegment(std::string &sequence, UserInputTeloscope userInput, uint64_t absPos);

    void insertWindowData(unsigned int seqPos, std::vector<WindowData>& pathWindows) {
        allWindows.push_back({seqPos, pathWindows});
    }

    void sortWindowsBySeqPos() {
        std::sort(allWindows.begin(), allWindows.end(), [](const std::pair<unsigned int, std::vector<WindowData>>& one, const std::pair<unsigned int, std::vector<WindowData>>& two) {
            return one.first < two.first;
        });
    }

    void printAllWindows() {
        
        std::cout << "Printing all windows finished!\n";

        // for (const auto& [seqPos, windows] : allWindows) {
        //     std::cout << "Sequence position: " << seqPos << "\n";

        //     for (const auto& window : windows) {
        //         std::cout << "Window start: " << window.windowStart << "\n";
        //         std::cout << "GC content: " << window.gcContent << "\n";
        //         std::cout << "Shannon entropy: " << window.shannonEntropy << "\n";

        //         for (const auto& [nucleotide, count] : window.nucleotideCounts) {
        //             std::cout << nucleotide << ": " << count << "\n";

        //         }
        //         for (const auto& [pattern, count] : window.patternCounts) {
        //             std::cout << pattern << ": " << count << "\n";
        //         }
        //     }
        // }
    }

    // void generateBEDFile(std::string outRoute) {
    //     std::string bedFileName = outRoute + "/teloscope.bed"; // Suffix and format change by data
    //     std::ofstream bedFile(bedFileName, std::ios::out);

    //     if (!bedFile.is_open()) {
    //         std::cerr << "Failed to open file: " << bedFileName << '\n';
    //         return;
    //     }

    //     for (const auto& [seqPos, windows] : allWindows) {
    //         for (const auto& window : windows) {
    //             for (const auto& [start, pattern] : window.patternBEDData) {
    //                 uint64_t end = start + pattern.size() - 1;
    //                 bedFile << cleanedHeader << "\t" << start << "\t" << end << "\t" << pattern << "\n";
    //             }
    //         }
    //     }

    //     bedFile.close();
    // }
};

#endif // TELOSCOPE_H/