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

    std::shared_ptr<TrieNode> getChild(const std::shared_ptr<TrieNode>& node, char ch) const {
        auto it = node->children.find(ch);
        if (it != node->children.end()) {
            return it->second;
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
    uint32_t blockDistance;
    uint16_t blockCounts;
};

struct WindowData {
    uint32_t windowStart;
    uint32_t currentWindowSize;
    uint32_t nucleotideCounts[4] = {0, 0, 0, 0};
    float gcContent;
    float shannonEntropy;
    // uint32_t winHDistance = 0;
    
    std::vector<uint32_t> winMatches;
    std::vector<uint8_t> hDistances; 
    uint16_t canonicalCounts = 0;
    uint16_t nonCanonicalCounts = 0;
    float canonicalDensity = 0.0f;
    float nonCanonicalDensity = 0.0f;
    bool hasCanDimer = false; 
    
    WindowData() : windowStart(0), currentWindowSize(0), gcContent(0.0f), shannonEntropy(0.0f), hasCanDimer(false) {}
};

struct SegmentData {
    std::vector<WindowData> windows;
    // std::vector<TelomereBlock> terminalBlocks;
    // std::vector<TelomereBlock> interstitialBlocks;
    std::unordered_map<std::string, std::vector<TelomereBlock>> mergedBlocks;
    std::vector<uint32_t> canonicalMatches;
    std::vector<uint32_t> nonCanonicalMatches;
    std::vector<uint32_t> segMatches;
    std::vector<TelomereBlock> segBlocks;
};


struct PathData {
    unsigned int seqPos;
    std::string header;
    std::vector<WindowData> windows; // Empty unless specified by user
    std::unordered_map<std::string, std::vector<TelomereBlock>> mergedBlocks;
    std::vector<uint32_t> canonicalMatches;
    std::vector<uint32_t> nonCanonicalMatches;
};


struct Stats {
    float min;
    float max;
    float mean;
    float median;
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


    Stats getStats(std::vector<float>& values);

public:

    Teloscope(UserInputTeloscope userInput) : userInput(userInput) {
        for (const auto& pattern : userInput.patterns) {
            trie.insertPattern(pattern);
        }
    }

    bool walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);

    void analyzeWindow(const std::string &window, uint32_t windowStart,
                        WindowData& windowData, WindowData& nextOverlapData,
                        SegmentData& segmentData, uint32_t segmentSize);

    SegmentData analyzeSegment(std::string &sequence, UserInputTeloscope userInput, uint64_t absPos);

    void insertWindowData(unsigned int seqPos, const std::string& header, std::vector<WindowData>& pathWindows);

    void sortBySeqPos();

    std::vector<TelomereBlock> getTelomereBlocks(const std::vector<uint32_t>& inputMatches, uint16_t mergeDist);

    std::vector<TelomereBlock> filterBlocks(const std::vector<TelomereBlock>& blocks);

    void writeBEDFile(std::ofstream& windowMetricsFile, std::ofstream& windowRepeatsFile,
                                std::ofstream& canonicalMatchFile, std::ofstream& noncanonicalMatchFile,
                                std::ofstream& allBlocksFile, std::ofstream& canonicalBlocksFile);

    void handleBEDFile();

    void printSummary();
};

#endif // TELOSCOPE_H/