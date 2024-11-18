#include "tools.h"

std::unordered_map<char, std::vector<char>> IUPAC = {
    {'A', {'A'}},
    {'C', {'C'}},
    {'G', {'G'}},
    {'T', {'T'}},
    {'R', {'A', 'G'}},
    {'Y', {'C', 'T'}},
    {'M', {'A', 'C'}},
    {'K', {'G', 'T'}},
    {'S', {'C', 'G'}},
    {'W', {'A', 'T'}},
    {'H', {'A', 'C', 'T'}},
    {'B', {'C', 'G', 'T'}},
    {'V', {'A', 'C', 'G'}},
    {'D', {'A', 'G', 'T'}},
    {'N', {'A', 'C', 'G', 'T'}}
};


void getCombinations(const std::string &pattern, std::string &current, size_t index, std::vector<std::string> &combinations) {
    if (index == pattern.size()) {
        combinations.push_back(current);
        return;
    }

    char base = pattern[index];
    for (char c : IUPAC[base]) {
        current[index] = c;
        getCombinations(pattern, current, index + 1, combinations);
    }
}


std::unordered_map<std::string, uint8_t> getHammingDistances(
    const std::vector<std::string> &patterns,
    const std::pair<std::string, std::string> &canonicalPatterns
) {
    // Helper lambda to get hDist
    auto hDist = [](const std::string &str1, const std::string &str2) -> uint8_t {
        if (str1.length() != str2.length()) {
            throw std::invalid_argument("Strings must be of equal length to compute Hamming distance.");
        }

        uint8_t distance = 0;
        for (size_t i = 0; i < str1.length(); ++i) {
            if (str1[i] != str2[i]) {
                ++distance;
            }
        }
        return distance;
    };

    std::unordered_map<std::string, uint8_t> hammingDistances;

    for (const auto &pattern : patterns) {
        uint8_t minDistance = std::min(hDist(pattern, canonicalPatterns.first), hDist(pattern, canonicalPatterns.second));
        hammingDistances[pattern] = minDistance;
    }

    return hammingDistances;
}