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
    uint32_t windowStart;
    uint32_t currentWindowSize;
    float gcContent;
    float shannonEntropy;
    std::unordered_map<char, uint64_t> nucleotideCounts;
    std::unordered_map<std::string, PatternData> patternMap; // Condensed pattern data
    
    WindowData() : windowStart(0), gcContent(0.0f), shannonEntropy(0.0f), nucleotideCounts{{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}} {}
};


class Teloscope {

    Trie trie; // Declare trie instance
    UserInputTeloscope userInput; // Declare user input instance
    std::vector<std::tuple<unsigned int, std::string, std::vector<WindowData>>> allWindows; // Assembly windows


    float getShannonEntropy(const std::unordered_map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);
    float getGCContent(const std::unordered_map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);


public:

    Teloscope(UserInputTeloscope userInput) : userInput(userInput) {
        for (const auto& pattern : userInput.patterns) {
        trie.insertPattern(pattern);
        }
    }

    bool walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps);

    void analyzeWindow(const std::string &window, uint32_t windowStart, WindowData& windowData);
    
    std::vector<WindowData> analyzeSegment(std::string &sequence, UserInputTeloscope userInput, uint64_t absPos);


    void insertWindowData(unsigned int seqPos, const std::string& header, std::vector<WindowData>& pathWindows) {
        allWindows.push_back(std::make_tuple(seqPos, header, pathWindows)); // Giulio: cleaner with struct
    }


    void sortWindowsBySeqPos() {
        std::sort(allWindows.begin(), allWindows.end(), [](const auto& one, const auto& two) {
            return std::get<0>(one) < std::get<0>(two);
        });
    }

    void printAllWindows() {
        std::cout << "Printing all windows in BEDs!\n";
    }


    // void generateBEDFile() {
    //     // Open GC/entropy files
    //     std::ofstream shannonFile(outRoute + "/shannonEntropy.bedgraph");
    //     std::ofstream gcContentFile(outRoute + "/gcContent.bedgraph");

    //     // Hold file streams for pattern data
    //     std::unordered_map<std::string, std::ofstream> patternMatchFiles;
    //     std::unordered_map<std::string, std::ofstream> patternCountFiles;
    //     std::unordered_map<std::string, std::ofstream> patternDensityFiles;

    //     // Open pattern-related files
    //     for (const auto& pattern : userInput.patterns) {
    //         patternMatchFiles[pattern].open(outRoute + "/" + pattern + "_matches.bed");
    //         patternCountFiles[pattern].open(outRoute + "/" + pattern + "_count.bedgraph");
    //         patternDensityFiles[pattern].open(outRoute + "/" + pattern + "_density.bedgraph");
    //     }

    //     // Write data for each window
    //     for (const auto& windowData : allWindows) {
    //         unsigned int seqPos;
    //         std::string header;
    //         std::vector<WindowData> windows;
    //         std::tie(seqPos, header, windows) = windowData; // Unpack the tuple

    //         for (const auto& window : windows) {
    //             uint32_t windowEnd = window.windowStart + window.currentWindowSize - 1;

    //             // Write window Shannon entropy and GC
    //             shannonFile << header << "\t" << window.windowStart << "\t"
    //                         << windowEnd << "\t"
    //                         << window.shannonEntropy << "\n";
    //             gcContentFile << header << "\t" << window.windowStart << "\t"
    //                         << windowEnd << "\t"
    //                         << window.gcContent << "\n";

    //             // Write pattern data
    //             for (const auto& [pattern, data] : window.patternMap) {
    //                 for (auto pos : data.positions) {
    //                     patternMatchFiles[pattern] << header << "\t"
    //                                             << window.windowStart + pos << "\t"
    //                                             << window.windowStart + pos + pattern.length() - 1 << "\t"
    //                                             << pattern << "\n";
    //                 }
    //                 patternCountFiles[pattern] << header << "\t" << window.windowStart << "\t"
    //                                         << windowEnd << "\t"
    //                                         << data.count << "\n";
    //                 patternDensityFiles[pattern] << header << "\t" << window.windowStart << "\t"
    //                                             << windowEnd << "\t"
    //                                             << data.density << "\n";
    //             }
    //         }
    //     }

    //     // Close all files
    //     shannonFile.close();
    //     gcContentFile.close();
    //     for (auto& [pattern, file] : patternMatchFiles) {
    //         file.close();
    //     }
    //     for (auto& [pattern, file] : patternCountFiles) {
    //         file.close();
    //     }
    //     for (auto& [pattern, file] : patternDensityFiles) {
    //         file.close();
    //     }
    // }

    void generateBEDFile() {
        std::ofstream shannonFile; // Declare file streams
        std::ofstream gcContentFile;
        
        std::unordered_map<std::string, std::ofstream> patternMatchFiles; // Hold file streams for pattern data
        std::unordered_map<std::string, std::ofstream> patternCountFiles;
        std::unordered_map<std::string, std::ofstream> patternDensityFiles;

        // Only create and write to files if their modes are enabled
        if (userInput.modeEntropy) {
            shannonFile.open(outRoute + "/shannonEntropy.bedgraph");
        }

        if (userInput.modeGC) {
            gcContentFile.open(outRoute + "/gcContent.bedgraph");
        }

        if (userInput.modeMatch) {
            for (const auto& pattern : userInput.patterns) {
                patternMatchFiles[pattern].open(outRoute + "/" + pattern + "_matches.bed");
                patternCountFiles[pattern].open(outRoute + "/" + pattern + "_count.bedgraph");
                patternDensityFiles[pattern].open(outRoute + "/" + pattern + "_density.bedgraph");
            }
        }

        // Write data for each window
        for (const auto& windowData : allWindows) {
            unsigned int seqPos;
            std::string header;
            std::vector<WindowData> windows;
            std::tie(seqPos, header, windows) = windowData; // Unpack the tuple

            for (const auto& window : windows) {
                uint32_t windowEnd = window.windowStart + window.currentWindowSize - 1;

                // Write window Shannon entropy if enabled
                if (userInput.modeEntropy) {
                    shannonFile << header << "\t" << window.windowStart << "\t"
                                << windowEnd << "\t"
                                << window.shannonEntropy << "\n";
                }

                // Write window GC content if enabled
                if (userInput.modeGC) {
                    gcContentFile << header << "\t" << window.windowStart << "\t"
                                << windowEnd << "\t"
                                << window.gcContent << "\n";
                }

                // Write pattern data if enabled
                if (userInput.modeMatch) {
                    for (const auto& [pattern, data] : window.patternMap) {
                        for (auto pos : data.positions) {
                            patternMatchFiles[pattern] << header << "\t"
                                                    << window.windowStart + pos << "\t"
                                                    << window.windowStart + pos + pattern.length() - 1 << "\t"
                                                    << pattern << "\n";
                        }
                        patternCountFiles[pattern] << header << "\t" << window.windowStart << "\t"
                                                << windowEnd << "\t"
                                                << data.count << "\n";
                        patternDensityFiles[pattern] << header << "\t" << window.windowStart << "\t"
                                                    << windowEnd << "\t"
                                                    << data.density << "\n";
                    }
                }
            }
        }

        // Close all files that were opened
        if (userInput.modeEntropy) {
            shannonFile.close();
        }
        if (userInput.modeGC) {
            gcContentFile.close();
        }
        if (userInput.modeMatch) {
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

    }
};

#endif // TELOSCOPE_H/