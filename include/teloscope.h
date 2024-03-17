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

void insertPattern(std::shared_ptr<TrieNode> root, const std::string &pattern);

void findPatternsInWindow(std::shared_ptr<TrieNode> root, const std::string &window, uint64_t windowStart, 
                        std::vector<std::tuple<uint64_t, std::string>> &patternBEDData, const UserInputTeloscope& userInput,
                        std::map<char, uint64_t> &nucleotideCounts, std::unordered_map<std::string, uint32_t> &patternCounts);

float getShannonEntropy(const std::map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);

float getGCContent(const std::map<char, uint64_t>& nucleotideCounts, uint32_t windowSize);

// template <typename T>
// void generateBEDFile(const std::string& header, const std::vector<std::tuple<uint64_t, T>>& data, 
//                     const std::string& fileName, const UserInputTeloscope& userInput, std::string &sequence, unsigned int pathId);

template <typename T>
void generateBEDFile(const std::string& header, const std::vector<std::tuple<uint64_t, T>>& data, 
                    const std::string& fileName, const UserInputTeloscope& userInput, 
                    std::string &sequence, unsigned int pathId, const std::map<unsigned int, uint64_t>& pathAbsPos);

void findTelomeres(std::string header, std::string &sequence, UserInputTeloscope userInput);

#endif // TELOSCOPE_H
