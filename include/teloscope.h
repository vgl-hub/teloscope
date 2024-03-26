#ifndef TELOSCOPE_H
#define TELOSCOPE_H

#include "input.h" // not in Mac's code
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
};


struct WindowData {
    uint64_t windowStart;
    float gcContent;
    float shannonEntropy;
    std::unordered_map<char, uint64_t> nucleotideCounts;
    std::unordered_map<std::string, uint32_t> patternCounts;
    std::vector<std::tuple<uint64_t, std::string>> patternBEDData;
    std::map<std::string, std::vector<std::tuple<uint64_t, uint64_t>>> patternCountData;
    std::map<std::string, std::vector<std::tuple<uint64_t, float>>> patternDensityData;

    
    WindowData() : windowStart(0), gcContent(0.0f), shannonEntropy(0.0f) {
        nucleotideCounts = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
    }
};

class Teloscope {

    Trie trie; // Declare trie instance
    UserInputTeloscope userInput; // Declare user input instance
    std::vector<std::pair<unsigned int, std::vector<WindowData>>> allWindows; // Assembly windows

    float getShannonEntropy(const std::unordered_map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);
    float getGCContent(const std::unordered_map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);


public:

    Teloscope(UserInputTeloscope userInput) : userInput(userInput) {}

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

    for (const auto& [seqPos, windows] : allWindows) {

        std::cout << "Sequence position: " << seqPos << "\n";

        for (const auto& window : windows) {
            std::cout << "Window start: " << window.windowStart << "\n";
            std::cout << "GC content: " << window.gcContent << "\n";
            std::cout << "Shannon entropy: " << window.shannonEntropy << "\n";

            for (const auto& [nucleotide, count] : window.nucleotideCounts) {
                std::cout << nucleotide << ": " << count << "\n";

            }
            for (const auto& [pattern, count] : window.patternCounts) {
                std::cout << pattern << ": " << count << "\n";
            }
        }
    }
}
};



// // Last stable
// void analyzeWindow(std::shared_ptr<TrieNode> root, const std::string &window, uint64_t windowStart, 
//                     std::vector<std::tuple<uint64_t, std::string>> &patternBEDData, const UserInputTeloscope& userInput,
//                     std::unordered_map<char, uint64_t> &nucleotideCounts, std::unordered_map<std::string, uint32_t> &patternCounts);


// void analyzeSegment(std::string header, std::string &sequence, UserInputTeloscope userInput, uint64_t absPos, unsigned int pathId);

// class BEDFileGenerator {
//     std::string header;
//     std::string fileName;
//     UserInputTeloscope userInput;
//     uint64_t absPos;

// public:
//     BEDFileGenerator(const std::string& header, const UserInputTeloscope& userInput, uint64_t absPos)
//         : header(header), userInput(userInput), absPos(absPos) {}

//     template <typename T>
//     void generateBEDFile(const std::vector<std::tuple<uint64_t, T>>& data, const std::string& fileName, const std::string& sequence);
// };


#endif // TELOSCOPE_H/
