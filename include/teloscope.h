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

    bool hasChild(const std::shared_ptr<TrieNode>& node, char ch) const { // merge with the following method
        return node->children.find(ch) != node->children.end();
    }

    std::shared_ptr<TrieNode> getChild(const std::shared_ptr<TrieNode>& node, char ch) const {
        if (hasChild(node, ch)) {
            return node->children[ch];
        }
        return nullptr; // Safe way to handle non-existent children
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
    std::vector<WindowData> allWindows;
    // std::vector<std::pair<unsigned int seqPos, std::vector<WindowData>>> allWindows;

    float getShannonEntropy(const std::unordered_map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);
    float getGCContent(const std::unordered_map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);


public:
    Teloscope(UserInputTeloscope userInput) : userInput(userInput) {}
    // Teloscope(const UserInputTeloscope& ui) : userInput(ui) {} // Initialize a Teloscope instance with user parameters
    // Teloscope teloscope(userInput);

    bool walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);
    void analyzeWindow(const std::string &window, uint64_t windowStart, WindowData& windowData);
    // std::vector<WindowData> analyzeSegment(std::string header, std::string &sequence, UserInputTeloscope userInput, uint64_t absPos, unsigned int pathId);
    std::vector<WindowData> analyzeSegment(std::string &sequence, UserInputTeloscope userInput, uint64_t absPos);

    void insertWindowData(std::vector<WindowData>& pathWindows) {
        allWindows.insert(allWindows.end(), pathWindows.begin(), pathWindows.end());
    }

    // void insertWindowData(unsigned int seqPos, std::vector<WindowData>& pathWindows) {
    //     allWindows.push_back({seqPos, pathWindows});
    // }

    // void sortWindowsBySeqPos() {
    //     std::sort(allWindows.begin(), allWindows.end(), [](const std::pair<unsigned int, std::vector<WindowData>>& one, const std::pair<unsigned int, std::vector<WindowData>>& two) {
    //         return one.first < two.first;
    //     });
    // }

    void printAllWindows() {
        for (const auto& window : allWindows) {
            std::cout << "Window start: " << window.windowStart << std::endl;
            std::cout << "GC content: " << window.gcContent << std::endl;
            std::cout << "Shannon entropy: " << window.shannonEntropy << std::endl;
            for (const auto& [nucleotide, count] : window.nucleotideCounts) {
                std::cout << nucleotide << ": " << count << std::endl;
            }
            for (const auto& [pattern, count] : window.patternCounts) {
                std::cout << pattern << ": " << count << std::endl;
            }
        }
    }

//     void printAllWindows() {
//     for (const auto& [seqPos, windows] : allWindows) {
//         std::cout << "Sequence position: " << seqPos << std::endl;
//         for (const auto& window : windows) {
//             std::cout << "Window start: " << window.windowStart << std::endl;
//             std::cout << "GC content: " << window.gcContent << std::endl;
//             std::cout << "Shannon entropy: " << window.shannonEntropy << std::endl;
//             for (const auto& [nucleotide, count] : window.nucleotideCounts) {
//                 std::cout << nucleotide << ": " << count << std::endl;
//             }
//             for (const auto& [pattern, count] : window.patternCounts) {
//                 std::cout << pattern << ": " << count << std::endl;
//             }
//         }
//     }
// }
};


// loop over a vector, for each element calls a function that prints the element

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
