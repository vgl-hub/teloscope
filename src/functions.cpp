#include "functions.h"
#include <unordered_map>

std::unordered_map<char, std::vector<char>> IUPAC_DNA_map = {
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


void generate_combinations(const std::string &pattern, std::string &current, size_t index, std::vector<std::string> &combinations) {
    if (index == pattern.size()) {
        combinations.push_back(current);
        return;
    }

    char base = pattern[index];
    for (char c : IUPAC_DNA_map[base]) {
        current[index] = c;
        generate_combinations(pattern, current, index + 1, combinations);
    }
}


