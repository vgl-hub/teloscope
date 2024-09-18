#ifndef TELOSCOPE_H
#define TELOSCOPE_H

#include "input.h" // not in Mac's code
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

// std::string removeCarriageReturns(const std::string& input);

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

    // std::shared_ptr<TrieNode> getChild(const std::shared_ptr<TrieNode>& node, char ch) const {
    //     auto it = node->children.find(ch);
    //     if (it != node->children.end()) {
    //         return it->second;
    //     }
    //     return nullptr;
    // }

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

    void annotateTelomeres(); // CHECK

    void writeBEDFile(std::ofstream& shannonFile, std::ofstream& gcContentFile,
                    // std::ofstream& telomereBEDFile, std::ofstream& telomereCountFile, // CHECK
                    std::unordered_map<std::string, std::ofstream>& patternMatchFiles,
                    std::unordered_map<std::string, std::ofstream>& patternCountFiles,
                    std::unordered_map<std::string, std::ofstream>& patternDensityFiles) {
        
        if (!userInput.keepWindowData) { // If windowData is not stored, return
            return;
        }

        for (const auto& windowData : allWindows) {
            unsigned int seqPos;
            std::string header;
            std::vector<WindowData> windows;
            std::tie(seqPos, header, windows) = windowData; // Unpack the tuple

            for (const auto& window : windows) {
                totalNWindows++; // Update total window count
                uint32_t windowEnd = window.windowStart + window.currentWindowSize; // Start is already 0-based

                // Write window Shannon entropy if enabled
                if (userInput.modeEntropy) {
                    shannonFile << header << "\t" << window.windowStart << "\t"
                                << windowEnd << "\t"
                                << window.shannonEntropy << "\n";
                    entropyValues.push_back(window.shannonEntropy); // Update entropy values
                }

                // Write window GC content if enabled
                if (userInput.modeGC) {
                    gcContentFile << header << "\t" << window.windowStart << "\t"
                                << windowEnd << "\t"
                                << window.gcContent << "\n";
                    gcContentValues.push_back(window.gcContent);
                }

                // Write pattern data if enabled
                if (userInput.modeMatch) {
                    for (const auto& [pattern, data] : window.patternMap) {
                        for (auto pos : data.wMatches) {
                            patternMatchFiles[pattern] << header << "\t"
                                                    << window.windowStart + pos << "\t"
                                                    << window.windowStart + pos + pattern.length() << "\t" // Start is already 0-based
                                                    << pattern << "\n";
                            patternCounts[pattern]++; // Update total pattern counts
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

            annotateTelomeres(); // CHECK

        }
    }


    void handleBEDFile() {
        std::ofstream shannonFile;
        std::ofstream gcContentFile;
        std::ofstream telomereBEDFile; // CHECK
        std::ofstream telomereCountFile; // CHECK
        std::unordered_map<std::string, std::ofstream> patternMatchFiles; // Jack: replace with vector to reduce cache locality?
        std::unordered_map<std::string, std::ofstream> patternCountFiles;
        std::unordered_map<std::string, std::ofstream> patternDensityFiles;
        std::cout << "Reporting window matches and metrics in BED/BEDgraphs...\n";

        // Open files once if their modes are enabled
        if (userInput.modeEntropy) {
            shannonFile.open(outRoute + "/shannonEntropy.bedgraph");
        }

        if (userInput.modeGC) {
            gcContentFile.open(outRoute + "/gcContent.bedgraph");
        }

        if (userInput.modeMatch) {
            // telomereBEDFile.open(outRoute + "/telomere_blocks.bed"); // CHECK
            // telomereCountFile.open(outRoute + "/telomere_block_counts.txt"); // CHECK

            for (const auto& pattern : userInput.patterns) {
                patternMatchFiles[pattern].open(outRoute + "/" + pattern + "_matches.bed");
                patternCountFiles[pattern].open(outRoute + "/" + pattern + "_count.bedgraph");
                patternDensityFiles[pattern].open(outRoute + "/" + pattern + "_density.bedgraph");
            }
        }

        // Write data for each window
        writeBEDFile(shannonFile, gcContentFile, 
                    // telomereBEDFile, telomereCountFile, // CHECK
                    patternMatchFiles, patternCountFiles, patternDensityFiles);

        // Close all files once
        if (userInput.modeEntropy) {
            shannonFile.close();
        }
        if (userInput.modeGC) {
            gcContentFile.close();
        }
        if (userInput.modeMatch) {
            // telomereBEDFile.close(); // CHECK
            // telomereCountFile.close(); // CHECK

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


    void printSummary() {
        std::cout << "\n+++Summary Report+++\n";
        std::cout << "Total windows analyzed:\t" << totalNWindows << "\n";
        std::cout << "Total input patterns found:\n";
        for (const auto& [pattern, count] : patternCounts) {
            std::cout << "Pattern:\t" << pattern << "\t" << count << "\n";
        }

        // For each pattern, print the path header with the highest number of matches - PENDING
        // For each pattern, print the path header with the lowest number of matches - PENDING
        if (userInput.keepWindowData) {
            std::cout << "Max Shannon Entropy:\t" << getMax(entropyValues) << "\n";
            std::cout << "Mean Shannon Entropy:\t" << getMean(entropyValues) << "\n";
            std::cout << "Median Shannon Entropy:\t" << getMedian(entropyValues) << "\n";
            std::cout << "Min Shannon Entropy:\t" << getMin(entropyValues) << "\n";

            std::cout << "Max GC Content:\t" << getMax(gcContentValues) << "\n";
            std::cout << "Mean GC Content:\t" << getMean(gcContentValues) << "\n";
            std::cout << "Median GC Content:\t" << getMedian(gcContentValues) << "\n";
            std::cout << "Min GC Content:\t" << getMin(gcContentValues) << "\n";
        }
    }
};

#endif // TELOSCOPE_H/