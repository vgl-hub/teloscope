#ifndef TOOLS_H
#define TOOLS_H

#include <stdint.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <stdexcept>

enum class ScaffoldType : uint8_t {
    T2T, GAPPED_T2T,
    MISASSEMBLY, GAPPED_MISASSEMBLY,
    INCOMPLETE, GAPPED_INCOMPLETE,
    NONE, GAPPED_NONE,
    DISCORDANT, GAPPED_DISCORDANT
};

const char* scaffoldTypeToString(ScaffoldType type);

struct Stats {
    float min = 0.0f;
    float max = 0.0f;
    float mean = 0.0f;
    float median = 0.0f;
};

Stats getStats(std::vector<float>& values);

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
