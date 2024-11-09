#ifndef TOOLS_H
#define TOOLS_H

#include <stdint.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <stdexcept>

void getCombinations(const std::string &pattern, std::string &current, size_t index, std::vector<std::string> &combinations);

std::unordered_map<std::string, uint8_t> getHammingDistances(
    const std::vector<std::string> &patterns,
    const std::pair<std::string, std::string> &canonicalPatterns
);

#endif // TOOLS_H
