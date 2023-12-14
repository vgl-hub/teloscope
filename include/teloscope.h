#include "input.h"

std::vector<uint64_t> getPatternFrequency(const std::vector<bool>& patternMatches, uint32_t windowSize, uint32_t step);

float getShannonEntropy(const std::string& window);

float getGCContent(const std::string& window);

void findTelomeres(std::string header, std::string &sequence, UserInputTeloscope userInput);