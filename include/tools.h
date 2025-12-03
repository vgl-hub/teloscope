#ifndef TOOLS_H
#define TOOLS_H

#include <stdint.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <stdexcept>
// #include <numeric>

void getCombinations(const std::string &pattern, std::string &current, size_t index, std::vector<std::string> &combinations);

std::vector<std::string> getEditVariants(const std::string &pattern, uint8_t maxDist);

std::unordered_map<std::string, uint8_t> getHammingDistances(
    const std::vector<std::string> &patterns,
    const std::pair<std::string, std::string> &canonicalPatterns
);

std::vector<std::string> expandPatterns(
    const std::vector<std::string> &rawPatterns,
    uint8_t editDistance);

std::vector<std::pair<std::string, bool>> expandPatternsWithOrientation(
    const std::vector<std::string> &rawPatterns,
    uint8_t editDistance,
    const std::string &canonicalFwd);

#endif // TOOLS_H
