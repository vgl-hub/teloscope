#ifndef TELOSCOPE_H
#define TELOSCOPE_H

#include "input.h"
#include "tools.h"
#include <iostream>
#include <map>
#include <stdint.h>
#include <vector>
#include <string>
#include <string_view>
#include <array>

class Trie {
    struct TrieNode {
        std::array<int32_t, 4> children = {-1, -1, -1, -1}; // A=0, C=1, G=2, T=3
        bool isEndOfWord = false;
        bool isForward = false; // true if pattern is forward-oriented (closer to canonicalFwd)
        bool isCanonical = false; // true if pattern is exactly canonicalFwd or canonicalRev
    };

    std::vector<TrieNode> nodes;          // contiguous node pool
    unsigned short int longestPatternSize = 0;

    // Map nucleotide to index (inline for speed)
    static int8_t charToIndex(char c) {
        switch (c) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default:  return -1; // skip non-ACGT
        }
    }

public:
    Trie() { nodes.emplace_back(); } // root at index 0

    void insertPattern(const std::string& pattern, bool isForward, bool isCanonical);

    // Return root index (always 0)
    int32_t getRoot() const { return 0; }

    // Fast child lookup: returns child index or -1
    int32_t getChild(int32_t nodeIdx, char ch) const {
        int8_t idx = charToIndex(ch);
        return (idx >= 0) ? nodes[nodeIdx].children[idx] : -1;
    }

    // Check if node marks end of a pattern
    bool isEnd(int32_t nodeIdx) const {
        return nodes[nodeIdx].isEndOfWord;
    }

    // Check if pattern ending at this node is forward-oriented
    bool isForward(int32_t nodeIdx) const {
        return nodes[nodeIdx].isForward;
    }

    // Check if pattern ending at this node is canonical (exact fwd or rev)
    bool isCanonical(int32_t nodeIdx) const {
        return nodes[nodeIdx].isCanonical;
    }

    unsigned short int getLongestPatternSize() const {
        return longestPatternSize;
    }
};


struct MatchInfo {
    bool isCanonical = false;
    bool isForward = false;
    uint64_t position = 0;
    uint16_t matchSize = 0;
    std::string matchSeq;
};


struct TelomereBlock {
    uint64_t start = 0;
    uint32_t blockLen = 0; // End = start + blockLen
    uint16_t blockCounts = 0;
    uint16_t forwardCount = 0;
    uint16_t reverseCount = 0;
    uint16_t canonicalCount = 0;
    uint16_t nonCanonicalCount = 0;
    uint32_t totalCovered = 0;
    uint32_t fwdCovered = 0;
    uint32_t canCovered = 0;
    bool hasValidOr = true;
    bool isLongest = false;
    char blockLabel = '\0'; // 'p', 'q', 'b' (balanced)
};

struct WindowData {
    uint64_t windowStart;
    uint32_t currentWindowSize;
    uint32_t nucleotideCounts[4] = {0, 0, 0, 0};
    float gcContent;
    float shannonEntropy;
    
    uint16_t canonicalCounts = 0;
    uint16_t nonCanonicalCounts = 0;
    uint16_t fwdCounts = 0;
    uint16_t revCounts = 0;
    float canonicalDensity = 0.0f;
    float nonCanonicalDensity = 0.0f;
    float fwdDensity = 0.0f;
    float revDensity = 0.0f;
    bool hasCanDimer = false; 
    
    WindowData() : windowStart(0), currentWindowSize(0), gcContent(0.0f), shannonEntropy(0.0f) {}
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
    ScaffoldType scaffoldType = ScaffoldType::NONE;
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
    uint32_t totalMisassembly = 0;
    uint32_t totalGappedMisassembly = 0;
    uint32_t totalIncomplete = 0;
    uint32_t totalGappedIncomplete = 0;
    uint32_t totalNone = 0;
    uint32_t totalGappedNone = 0;
    uint32_t totalDiscordant = 0;
    uint32_t totalGappedDiscordant = 0;

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


    static char computeBlockLabel(uint16_t forwardCount, uint16_t blockCounts);

    void computeSummaryCounts();

public:

    Teloscope(UserInputTeloscope userInput) : userInput(userInput) {
        for (const auto& [pattern, isForward] : this->userInput.patternInfo) {
            bool isCanonical = (pattern == this->userInput.canonicalFwd ||
                                pattern == this->userInput.canonicalRev);
            trie.insertPattern(pattern, isForward, isCanonical);
        }
    }

    bool walkSegment(InSegment* segment, InSequences& inSequences);

    bool walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);

    void analyzeWindow(const std::string_view &window, uint64_t windowStart,
                        WindowData& windowData, WindowData& nextOverlapData,
                        SegmentData& segmentData, uint64_t segmentSize, uint64_t absPos);

    SegmentData scanSegment(std::string &sequence, uint64_t absPos, bool tipsOnly);

    void sortBySeqPos();

    std::vector<TelomereBlock> getTeloBlocks(
        std::vector<MatchInfo>& matches, 
        uint16_t mergeDist,
        bool needsSorting = false,
        std::vector<MatchInfo>* recycleTarget = nullptr,
        bool recycleToStart = false);

    std::vector<TelomereBlock> extendBlocks(std::vector<TelomereBlock> &blocks,
    uint16_t maxBlockDist, float densityCutoff, uint64_t segmentSize, uint64_t absPos);

    void labelTerminalBlocks(std::vector<TelomereBlock>& blocks, uint16_t gaps,
                        std::string& terminalLabel, ScaffoldType& scaffoldType,
                        uint64_t pathSize, uint32_t terminalLimit);
    
    std::vector<TelomereBlock> filterITSBlocks(const std::vector<TelomereBlock>& interstitialBlocks);
    
    std::string getChrType(const std::string& labels, uint16_t gaps);
    
    void writeBEDFile(std::ofstream& windowDensityFile,
                    std::ofstream& windowCanonicalRatioFile,
                    std::ofstream& windowStrandRatioFile,
                    std::ofstream& windowGCFile,
                    std::ofstream& windowEntropyFile,
                    std::ofstream& canonicalMatchFile,
                    std::ofstream& noncanonicalMatchFile,
                    std::ofstream& terminalBlocksFile,
                    std::ofstream& interstitialBlocksFile);


    void handleBEDFile();

    void printSummary();
};

#endif // TELOSCOPE_H/