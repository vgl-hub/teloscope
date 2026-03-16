#ifndef TELOSCOPE_H
#define TELOSCOPE_H

#include "input.h" // not in Mac's code
#include <iostream>
#include <map>
#include <stdint.h>
#include <vector>
#include <string>
#include <string_view>
#include <array>

enum class ScaffoldType : uint8_t {
    T2T, GAPPED_T2T,
    MISSASSEMBLY, GAPPED_MISSASSEMBLY,
    INCOMPLETE, GAPPED_INCOMPLETE,
    NONE, GAPPED_NONE,
    DISCORDANT, GAPPED_DISCORDANT
};

const char* scaffoldTypeToString(ScaffoldType type);

class Trie {
    struct TrieNode {
        std::array<int32_t, 4> children = {-1, -1, -1, -1}; // A=0, C=1, G=2, T=3
        bool isEndOfWord = false;
        bool isForward = false; // true if pattern is forward-oriented (closer to canonicalFwd)
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

    void insertPattern(const std::string& pattern, bool isForward);

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

    unsigned short int getLongestPatternSize() const {
        return longestPatternSize;
    }
};


struct MatchInfo {
    bool isCanonical = false;
    bool isForward = false;
    uint32_t position = 0;
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
    uint32_t windowStart;
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
    ScaffoldType scaffoldType = ScaffoldType::NONE;
};


struct Stats {
    float min = 0.0f;
    float max = 0.0f;
    float mean = 0.0f;
    float median = 0.0f;
};


class Teloscope {
    Trie trie; // Declare trie instance
    UserInputTeloscope userInput; // Declare user input instance
    std::vector<PathData> allPathData; // Assembly data
    
    // Cached string_views for canonical patterns (avoid reconstruction per window)
    std::string_view canonicalFwdView;
    std::string_view canonicalRevView;

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


    Stats getStats(std::vector<float>& values);

    static char computeBlockLabel(uint16_t forwardCount, uint16_t blockCounts);

    void computeSummaryCounts();

public:

    Teloscope(UserInputTeloscope userInput) : userInput(userInput),
        canonicalFwdView(this->userInput.canonicalFwd),
        canonicalRevView(this->userInput.canonicalRev) {
        for (const auto& [pattern, isForward] : this->userInput.patternInfo) {
            trie.insertPattern(pattern, isForward);
        }
    }

    bool walkSegment(InSegment* segment, InSequences& inSequences);

    bool walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);

    void analyzeWindow(const std::string_view &window, uint32_t windowStart,
                        WindowData& windowData, WindowData& nextOverlapData,
                        SegmentData& segmentData, uint32_t segmentSize, uint32_t absPos);

    SegmentData scanSegment(std::string &sequence, uint32_t absPos, bool tipsOnly);

    void sortBySeqPos();

    std::vector<TelomereBlock> getTeloBlocks(
        std::vector<MatchInfo>& matches, 
        uint16_t mergeDist,
        bool needsSorting = false,
        std::vector<MatchInfo>* recycleTarget = nullptr,
        bool recycleToStart = false);

    std::vector<TelomereBlock> extendBlocks(std::vector<TelomereBlock> &blocks, 
    uint16_t maxBlockDist, float densityCutoff, uint32_t segmentSize, uint32_t absPos);

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