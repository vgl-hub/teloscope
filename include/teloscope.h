#ifndef TELOSCOPE_H
#define TELOSCOPE_H

#include "input.h"
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

// Forward declaration
struct TrieNode;

float getShannonEntropy(const std::map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);
float getGCContent(const std::map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);
void insertPattern(std::shared_ptr<TrieNode> root, const std::string &pattern);
void findPatternsInWindow(std::shared_ptr<TrieNode> root, const std::string &window,
                        uint64_t windowStart, std::vector<std::tuple<uint64_t, std::string>> &patternBEDData,
                        std::map<std::string, uint64_t> &lastPatternPositions, uint32_t step,
                        std::map<char, uint64_t> &nucleotideCounts, std::map<char, uint64_t> &prevOverlapCounts);
template <typename T>
void generateBEDFile(const std::string& header, const std::vector<std::tuple<uint64_t, T>>& data, 
                    const std::string& fileName, uint32_t windowSize);
void findTelomeres(std::string header, std::string &sequence, UserInputTeloscope userInput);

#endif // TELOSCOPE_H
