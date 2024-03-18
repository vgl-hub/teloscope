#ifndef TELOSCOPE_H
#define TELOSCOPE_H

#include "input.h"
#include <map>
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

public:

    Trie() : root(std::make_shared<TrieNode>()) {}

    void insertPattern(const std::string& pattern);

    void findPatternsInWindow(const std::string &window, uint64_t windowStart,
                            std::vector<std::tuple<uint64_t, std::string>> &patternBEDData, const UserInputTeloscope& userInput,
                            std::map<char, uint64_t> &nucleotideCounts, std::unordered_map<std::string, uint32_t> &patternCounts);
};

// void findPatternsInWindow(std::shared_ptr<TrieNode> root, const std::string &window, uint64_t windowStart, 
//                         std::vector<std::tuple<uint64_t, std::string>> &patternBEDData, const UserInputTeloscope& userInput,
//                         std::map<char, uint64_t> &nucleotideCounts, std::unordered_map<std::string, uint32_t> &patternCounts);

float getShannonEntropy(const std::map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);

float getGCContent(const std::map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);

template <typename T>
void generateBEDFile(const std::string& header, const std::vector<std::tuple<uint64_t, T>>& data, 
                    const std::string& fileName, const UserInputTeloscope& userInput, std::string &sequence);

void findTelomeres(std::string header, std::string &sequence, UserInputTeloscope userInput);

#endif // TELOSCOPE_H
