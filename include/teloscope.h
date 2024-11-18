#ifndef TELOSCOPE_H
#define TELOSCOPE_H

#include "input.h" // not in Mac's code
#include <iostream>
#include <map>
#include <stdint.h>
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

    bool hasChild(const std::shared_ptr<TrieNode>& node, char ch) const { // Merge with the following method
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


struct TelomereBlock {
    uint64_t start;
    uint16_t blockLen; // End = start + blockLen
};

struct WindowData {
    uint32_t windowStart;
    uint32_t currentWindowSize;
    uint32_t nucleotideCounts[4] = {0, 0, 0, 0};
    float gcContent;
    float shannonEntropy;
    // uint32_t winHDistance = 0;
    
    std::vector<TelomereBlock> winBlocks;
    std::vector<uint8_t> hDistances; 

    std::vector<uint32_t> canonicalMatches;
    std::vector<uint32_t> nonCanonicalMatches;
    std::vector<uint32_t> windowMatches;
    uint16_t canonicalCounts = 0;
    uint16_t nonCanonicalCounts = 0;
    float canonicalDensity = 0.0f;
    float nonCanonicalDensity = 0.0f;
    
    WindowData() : windowStart(0), currentWindowSize(0), gcContent(0.0f), shannonEntropy(0.0f) {}
};

struct SegmentData {
    std::vector<WindowData> windows;
    std::unordered_map<std::string, std::vector<TelomereBlock>> mergedBlocks;
};


struct PathData {
    unsigned int seqPos;
    std::string header;
    std::vector<WindowData> windows; // Empty unless specified by user
    std::unordered_map<std::string, std::vector<TelomereBlock>> mergedBlocks;
};


class Teloscope {
    Trie trie; // Declare trie instance
    UserInputTeloscope userInput; // Declare user input instance
    std::vector<PathData> allPathData; // Assembly data

    int totalNWindows = 0; // Total windows analyzed
    std::unordered_map<std::string, int> patternCounts; // Total counts
    std::vector<float> entropyValues; // Total entropy values
    std::vector<float> gcContentValues; // Total GC content values


    inline float getShannonEntropy(const uint32_t nucleotideCounts[4], uint32_t windowSize) {
        float entropy = 0.0;
        for (int i = 0; i < 4; ++i) {
            if (nucleotideCounts[i] > 0) {
                float probability = static_cast<float>(nucleotideCounts[i]) / windowSize;
                entropy -= probability * std::log2(probability);
            }
        }
        return entropy;
    }


    inline float getGCContent(const uint32_t nucleotideCounts[4], uint32_t windowSize) {
        uint32_t gcCount = nucleotideCounts[1] + nucleotideCounts[2]; // Indices: 1 = C, 2 = G
        return static_cast<float>(gcCount) / windowSize * 100.0;
    }

    float getMean(const std::vector<float>& values);
    float getMedian(std::vector<float> values);
    float getMin(const std::vector<float> values);
    float getMax(const std::vector<float> values);

public:

    Teloscope(UserInputTeloscope userInput) : userInput(userInput) {
        for (const auto& pattern : userInput.patterns) {
            trie.insertPattern(pattern);
        }
    }

    bool walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);

    void analyzeWindow(const std::string &window, uint32_t windowStart, WindowData& windowData, WindowData& nextOverlapData); 

    SegmentData analyzeSegment(std::string &sequence, UserInputTeloscope userInput, uint64_t absPos);

    void insertWindowData(unsigned int seqPos, const std::string& header, std::vector<WindowData>& pathWindows);

    void sortBySeqPos();

    // std::vector<TelomereBlock> getTelomereBlocks(const std::vector<uint32_t>& inputMatches, uint64_t windowStart);
    std::vector<TelomereBlock> getTelomereBlocks(const std::vector<uint32_t>& inputMatches, uint64_t windowStart, uint32_t currentWindowSize);

    std::vector<TelomereBlock> mergeTelomereBlocks(const std::vector<TelomereBlock>& winBlocks);

    void writeBEDFile(std::ofstream& shannonFile, std::ofstream& gcContentFile,
                                std::ofstream& canonicalMatchFile, std::ofstream& noncanonicalMatchFile,
                                std::ofstream& windowRepeatsFile,
                                std::ofstream& allBlocksFile, std::ofstream& canonicalBlocksFile, std::ofstream& noncanonicalBlocksFile);


    void handleBEDFile();

    void printSummary();
};

#endif // TELOSCOPE_H/