#include "input.h"

float getShannonEntropy(const std::string& window);

float getGCContent(const std::string& window);

template <typename T>
void generateBEDFile(const std::string& header, const std::vector<std::tuple<uint64_t, T>>& data, const std::string& fileName, uint32_t windowSize);

void findTelomeres(std::string header, std::string &sequence, UserInputTeloscope userInput);