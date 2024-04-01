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
    //         return it->second; // Return shared_ptr to the TrieNode if found
    //     }
    //     return nullptr; // If character is not a child of the given node
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
        
        std::cout << "Printing all windows works!\n";

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

    void generateBEDFile() {
        // Open files
        std::ofstream shannonFile(outRoute + "/shannonEntropy.bedgraph");
        std::ofstream gcContentFile(outRoute + "/gcContent.bedgraph");

        // Hold file streams for pattern data
        std::unordered_map<std::string, std::ofstream> patternMatchFiles;
        std::unordered_map<std::string, std::ofstream> patternCountFiles;
        std::unordered_map<std::string, std::ofstream> patternDensityFiles;

        // Open files for all patterns
        for (const auto& pattern : userInput.patterns) {
            patternMatchFiles[pattern].open(outRoute + "/" + pattern + "_matches.bed");
            patternCountFiles[pattern].open(outRoute + "/" + pattern + "_count.bedgraph");
            patternDensityFiles[pattern].open(outRoute + "/" + pattern + "_density.bedgraph");
        }

        // Write data for each window
        for (const auto& [seqPos, windows] : allWindows) {
            for (const auto& window : windows) {

                // std::string cleanHeader = cleanString(header); // Jack: header = cleanString(path->getHeader());
                std::string cleanHeader = "header";

                uint64_t windowEnd = window.windowStart + userInput.windowSize - 1;

                // Write Shannon entropy and GC content
                shannonFile << cleanHeader << "\t" << window.windowStart << "\t"
                            << windowEnd << "\t"
                            << window.shannonEntropy << "\n";
                gcContentFile << cleanHeader << "\t" << window.windowStart << "\t"
                            << windowEnd << "\t"
                            << window.gcContent << "\n";

                // Write pattern data
                for (const auto& [pattern, data] : window.patternMap) {
                    for (auto pos : data.positions) {
                        patternMatchFiles[pattern] << cleanHeader << "\t"
                                                << window.windowStart + pos << "\t"
                                                << window.windowStart + pos + pattern.length() - 1 << "\t"
                                                << pattern << "\n";
                    }
                    patternCountFiles[pattern] << cleanHeader << "\t" << window.windowStart << "\t"
                                            << windowEnd << "\t"
                                            << data.count << "\n";
                    patternDensityFiles[pattern] << cleanHeader << "\t" << window.windowStart << "\t"
                                                << windowEnd << "\t"
                                                << data.density << "\n";
                }
            }
        }

        // Close all files
        shannonFile.close();
        gcContentFile.close();
        for (auto& [pattern, file] : patternMatchFiles) {
            file.close();
        }
        for (auto& [pattern, file] : patternCountFiles) {
            file.close();
        }
        for (auto& [pattern, file] : patternDensityFiles) {
            file.close();
        }
    }

};

#endif // TELOSCOPE_H/