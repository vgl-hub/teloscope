#include "tools.h"
#include "functions.h"

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


std::vector<std::string> getEditVariants(const std::string &pattern, uint8_t maxDist) {
    // Generate all substitution variants up to maxDist edit distance. 
    std::vector<std::string> variants;
    
    if (maxDist == 0) {
        return variants; // No variants needed for exact matching
    }
    
    const char nucleotides[] = {'A', 'C', 'G', 'T'};
    
    // For edit distance 1: generate all single-substitution variants
    if (maxDist >= 1) {
        for (size_t i = 0; i < pattern.size(); ++i) {
            char original = pattern[i];
            
            // "N" or other IUPAC codes skipped.
            if (original != 'A' && original != 'C' && original != 'G' && original != 'T') {
                continue;
            }
            
            for (char c : nucleotides) {
                if (c != original) {
                    std::string variant = pattern;
                    variant[i] = c;
                    variants.push_back(variant);
                }
            }
        }
    }
    
    // For edit distance 2+: generate from distance-1 variants (iterate by index to avoid copy)
    if (maxDist >= 2) {
        size_t dist1End = variants.size(); // Snapshot count before adding distance-2
        for (size_t v = 0; v < dist1End; ++v) {
            std::vector<std::string> d2vars = getEditVariants(variants[v], 1);
            variants.insert(variants.end(), d2vars.begin(), d2vars.end());
        }
    }
    
    return variants;
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


std::vector<std::string> expandPatterns(
    const std::vector<std::string> &rawPatterns,
    uint8_t editDistance) {

    std::vector<std::string> seeds = rawPatterns;

    std::vector<std::string> preparedPatterns;

    for (const auto &seed : seeds) {
        if (seed.empty()) {
            continue;
        }

        std::vector<std::string> combinations;
        std::string working = seed;
        getCombinations(seed, working, 0, combinations);

        for (const auto &combo : combinations) {
            std::vector<std::string> variants;
            variants.push_back(combo);

            if (editDistance > 0) {
                auto edits = getEditVariants(combo, editDistance);
                variants.insert(variants.end(), edits.begin(), edits.end());
            }

            for (const auto &variant : variants) {
                preparedPatterns.emplace_back(variant);
                preparedPatterns.emplace_back(revCom(variant));
            }
        }
    }

    std::sort(preparedPatterns.begin(), preparedPatterns.end());
    preparedPatterns.erase(
        std::unique(preparedPatterns.begin(), preparedPatterns.end()),
        preparedPatterns.end());

    return preparedPatterns;
}