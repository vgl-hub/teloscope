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

std::unordered_map<std::string, uint8_t> getHammingDistances(
    const std::vector<std::string> &patterns,
    const std::pair<std::string, std::string> &canonicalPatterns
);

// inline float getMean(const std::vector<float>& values) {
//     if (values.empty()) return 0.0;
//     float sum = std::accumulate(values.begin(), values.end(), 0.0);
//     return sum / values.size();
// }

// inline float getMedian(std::vector<float> values) {
//     if (values.empty()) return 0.0;
//     std::sort(values.begin(), values.end());
//     size_t size = values.size();
//     if (size % 2 == 0) {
//         return (values[size / 2 - 1] + values[size / 2]) / 2;
//     } else {
//         return values[size / 2];
//     }
// }

// inline float getMin(const std::vector<float>& values) {
//     if (values.empty()) return 0.0;
//     return *std::min_element(values.begin(), values.end());
// }

// inline float getMax(const std::vector<float>& values) {
//     if (values.empty()) return 0.0;
//     return *std::max_element(values.begin(), values.end());
// }


#endif // TOOLS_H
