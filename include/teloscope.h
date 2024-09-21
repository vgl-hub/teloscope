#ifndef TELOSCOPE_H
#define TELOSCOPE_H

#include "input.h" // not in Mac's code
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

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

    unsigned short int getLongestPatternSize() const {
        return longestPatternSize;
    }
};


struct PatternData {
    std::vector<uint32_t> wMatches; // Match indexes in window
    uint32_t count = 0; // Total pattern count
    float density = 0.0f; // Density of the pattern
};


struct WindowData {
    uint32_t windowStart;
    uint32_t currentWindowSize;
    float gcContent;
    float shannonEntropy;
    std::unordered_map<char, uint32_t> nucleotideCounts;
    std::unordered_map<std::string, PatternData> patternMap; // Condensed pattern data
    
    WindowData() : windowStart(0), gcContent(0.0f), shannonEntropy(0.0f), nucleotideCounts{{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}} {}
};

struct TelomereBlock {
    uint64_t start;
    uint64_t end;
};

class Teloscope {

    Trie trie; // Declare trie instance
    UserInputTeloscope userInput; // Declare user input instance
    std::vector<std::tuple<unsigned int, std::string, std::vector<WindowData>>> allWindows; // Assembly windows
    std::vector<std::tuple<unsigned int, std::string, std::vector<TelomereBlock>>> allTelomereBlocks; // Assembly elomere blocks

    int totalNWindows = 0; // Total windows analyzed
    std::unordered_map<std::string, int> patternCounts; // Total counts
    std::vector<float> entropyValues; // Total entropy values
    std::vector<float> gcContentValues; // Total GC content values

    float getShannonEntropy(const std::unordered_map<char, uint32_t>& nucleotideCounts, uint32_t windowSize);
    float getGCContent(const std::unordered_map<char, uint32_t>& nucleotideCounts, uint32_t windowSize);
    void getPatternDensities(WindowData& windowData, uint32_t windowSize);

    float getMean(const std::vector<float>& values);
    float getMedian(std::vector<float> values);
    float getMin(std::vector<float> values);
    float getMax(std::vector<float> values);

public:

    Teloscope(UserInputTeloscope userInput) : userInput(userInput) {
        for (const auto& pattern : userInput.patterns) {
            trie.insertPattern(pattern);
        }
    }

    bool walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);

    void analyzeWindow(const std::string &window, uint32_t windowStart, WindowData& windowData);    

    std::vector<WindowData> analyzeSegment(std::string &sequence, UserInputTeloscope userInput, uint64_t absPos);

    void insertWindowData(unsigned int seqPos, const std::string& header, std::vector<WindowData>& pathWindows);

    void sortWindowsBySeqPos();

    void writeBEDFile(std::ofstream& shannonFile, std::ofstream& gcContentFile,
                    std::unordered_map<std::string, std::ofstream>& patternMatchFiles,
                    std::unordered_map<std::string, std::ofstream>& patternCountFiles,
                    std::unordered_map<std::string, std::ofstream>& patternDensityFiles);

    void handleBEDFile();

    void printSummary();
};

#endif // TELOSCOPE_H/