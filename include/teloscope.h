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


struct MatchInfo {
    bool isCanonical;
    bool isForward;
    uint32_t position;
    uint16_t matchSize;
    std::string matchSeq;

    MatchInfo() : isCanonical(false), isForward(false), matchSize(0) {}
};


struct TelomereBlock {
    uint64_t start;
    uint32_t blockLen; // End = start + blockLen
    uint16_t blockCounts;
    uint16_t forwardCount;
    uint16_t reverseCount;
    uint16_t canonicalCount;
    uint16_t nonCanonicalCount;
    uint32_t totalCovered;
    uint32_t fwdCovered;
    uint32_t canCovered;
    bool hasValidOr;
    bool isLongest;
    char blockLabel; // 'p', 'q', 'u'

    TelomereBlock() : hasValidOr(true), isLongest(false) {}
};

struct WindowData {
    uint32_t windowStart;
    uint32_t currentWindowSize;
    uint32_t nucleotideCounts[4] = {0, 0, 0, 0};
    float gcContent;
    float shannonEntropy;
    // uint32_t winHDistance = 0;
    
    std::vector<MatchInfo> terminalFwdMatches;
    std::vector<MatchInfo> terminalRevMatches;
    std::vector<MatchInfo> interstitialMatches;
    std::vector<uint8_t> hDistances; 
    uint16_t canonicalCounts = 0;
    uint16_t nonCanonicalCounts = 0;
    uint16_t fwdCounts = 0;
    uint16_t revCounts = 0;
    float canonicalDensity = 0.0f;
    float nonCanonicalDensity = 0.0f;
    float fwdDensity = 0.0f;
    float revDensity = 0.0f;
    bool hasCanDimer = false; 
    
    WindowData() : windowStart(0), currentWindowSize(0), gcContent(0.0f), shannonEntropy(0.0f), hasCanDimer(false) {}
};

struct SegmentData {
    std::vector<WindowData> windows;
    std::vector<TelomereBlock> terminalBlocks;
    std::vector<TelomereBlock> interstitialBlocks;
    std::vector<MatchInfo> canonicalMatches;
    std::vector<MatchInfo> nonCanonicalMatches;
    std::vector<MatchInfo> terminalFwdMatches;
    std::vector<MatchInfo> terminalRevMatches;
    std::vector<MatchInfo> interstitialMatches;
};


struct PathData {
    unsigned int seqPos;
    std::string header;
    uint16_t gaps = 0;
    uint64_t pathSize;
    std::vector<WindowData> windows;
    std::vector<TelomereBlock> terminalBlocks;
    std::vector<TelomereBlock> interstitialBlocks;
    std::vector<MatchInfo> canonicalMatches;
    std::vector<MatchInfo> nonCanonicalMatches;
    std::string terminalLabel;
    std::string scaffoldType;
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

    // Assembly Summary 
    uint32_t totalPaths = 0;
    uint32_t totalNWindows = 0;
    uint32_t totalTelomeres = 0;
    uint32_t totalITS = 0;
    uint32_t totalCanMatches = 0;
    uint32_t totalGaps = 0;

    // Telomere stats
    float teloMean = 0.0f;
    float teloMedian = 0.0f;
    float teloMin = 0.0f;
    float teloMax = 0.0f;

    // Chr/scaffold type summary
    uint32_t totalT2T = 0;
    uint32_t totalGappedT2T = 0;
    uint32_t totalMissassembly = 0;
    uint32_t totalGappedMissassembly = 0;
    uint32_t totalIncomplete = 0;
    uint32_t totalGappedIncomplete = 0;
    uint32_t totalNone = 0;
    uint32_t totalGappedNone = 0;
    uint32_t totalErrors = 0;
    uint32_t totalGappedErrors = 0;

    inline float getShannonEntropy(const uint32_t nucleotideCounts[4], uint32_t windowSize) {
        float entropy = 0.0;
        for (int i = 0; i < 4; ++i) {
            if (nucleotideCounts[i] > 0) {
                float probability = static_cast<float>(nucleotideCounts[i]) / windowSize;
                entropy -= probability * std::log2(probability);
            }
        }
        return std::round(entropy * 1000.0f) / 1000.0f; // Round to 3 decimal places
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

    bool walkSegment(InSegment* segment, InSequences& inSequences);

    bool walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);

    void analyzeWindow(const std::string_view &window, uint32_t windowStart,
                        WindowData& windowData, WindowData& nextOverlapData,
                        SegmentData& segmentData, uint32_t segmentSize, uint32_t absPos);

    SegmentData analyzeSegment(std::string &sequence, UserInputTeloscope userInput, uint32_t absPos);

    SegmentData analyzeSegmentTips(std::string &sequence, UserInputTeloscope &userInput, uint32_t absPos);

    void sortBySeqPos();

    std::vector<TelomereBlock> getBlocksRecycle(
        const std::vector<MatchInfo>& matches, 
        uint16_t mergeDist,
        std::vector<MatchInfo>& interstitialMatches,
        bool recycleToStart);

        std::vector<TelomereBlock> getBlocks(
            std::vector<MatchInfo>& matches, 
            uint16_t mergeDist, bool needsSorting);

    std::vector<TelomereBlock> extendBlocks(std::vector<TelomereBlock> &blocks, 
    uint16_t maxBlockDist, float densityCutoff, uint32_t segmentSize, uint32_t absPos);

    void labelTerminalBlocks(std::vector<TelomereBlock>& blocks, uint16_t gaps,
                        std::string& terminalLabel, std::string& scaffoldType);

    // std::vector<TelomereBlock> filterTerminalBlocks(const std::vector<TelomereBlock>& blocks);
    
    std::vector<TelomereBlock> filterITSBlocks(const std::vector<TelomereBlock>& interstitialBlocks);
    
    std::string getChrType(const std::string& labels, uint16_t gaps);
    
    void writeBEDFile(std::ofstream& windowMetricsFile,
                    std::ofstream& canonicalMatchFile,
                    std::ofstream& noncanonicalMatchFile,
                    std::ofstream& terminalBlocksFile,
                    std::ofstream& interstitialBlocksFile);


    void handleBEDFile();

    void printSummary();
};

#endif // TELOSCOPE_H/