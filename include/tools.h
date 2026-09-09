#ifndef TOOLS_H
#define TOOLS_H

#include <stdint.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <stdexcept>

// completeness only: how many arms. Gappedness comes from the gap list, and whether
// the arms are plausible is a separate axis, so neither can consume this answer.
enum class ScaffoldType : uint8_t {
    T2T, INCOMPLETE, NONE
};

// plausibility, accumulated as a set so a scaffold can carry more than one
enum AnomalyFlag : uint8_t {
    ANOM_DISC_P = 1,
    ANOM_DISC_Q = 2,
    ANOM_BAL_P  = 4,
    ANOM_BAL_Q  = 8,
    ANOM_EXTRA  = 16
};

const char* scaffoldTypeToString(ScaffoldType type);
std::string anomalyFlagsToString(uint8_t flags);

struct Stats {
    float min = 0.0f;
    float max = 0.0f;
    float mean = 0.0f;
    float median = 0.0f;
};

Stats getStats(std::vector<float>& values);

constexpr size_t maxCombinations = 4096;

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
