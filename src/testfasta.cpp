#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <algorithm>

std::vector<std::string> askForPatterns() {
    std::vector<std::string> patterns;
    std::string input;
    std::cout << "Enter telomeric patterns (type 'done' when finished, default is 'TTAGGG' if none provided): ";
    while (true) {
        std::cin >> input;
        if (input == "done") {
            break;
        }
        patterns.push_back(input);
    }

    if (patterns.empty()) {
        patterns.push_back("TTAGGG");
    }
    return patterns;
}
// std::pair
std::vector<int> askForRepeats(int numPatterns) {
    std::vector<int> repeats;
    int input;
    std::cout << "Enter the number of repeats for each pattern (type '-1' when finished): ";
    for (int i = 0; i < numPatterns; ++i) {
        while (true) {
            std::cin >> input;
            if (input == -1) {
                break;
            }
            if (input > 0) {
                repeats.push_back(input);
                break;
            }
            std::cout << "Please enter a positive integer." << std::endl;
        }
        if (input == -1) break;
    }
    return repeats;
}

char askForMode() {
    char mode;
    while (true) {
        std::cout << "Choose mode (c for concatenate, r for random): ";
        std::cin >> mode;
        if (mode == 'c' || mode == 'r') {
            break;
        }
        std::cout << "Invalid input. Please enter 'c' or 'r'." << std::endl;
    }
    return mode;
}
// unlocalized vs unplaced
int askForMismatches() {
    int mismatches;
    while (true) {
        std::cout << "Enter number of mismatches: ";
        if (std::cin >> mismatches && mismatches >= 0) {
            break;
        }
        std::cout << "Please enter a non-negative integer." << std::endl;
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }  
    return mismatches;
}

std::string generateSequence(const std::vector<std::string>& patterns, 
                             const std::vector<int>& repeats, 
                             char mode, int mismatches) {
    std::string sequence = "";
    // Concatenate mode
    if (mode == 'c') {
        for (size_t i = 0; i < patterns.size(); ++i) {
            for (int j = 0; j < repeats[i]; ++j) {
                sequence += patterns[i];
            }
        }
    } else if (mode == 'r') { // Random mode
        std::vector<std::string> blocks;
        for (size_t i = 0; i < patterns.size(); ++i) {
            for (int j = 0; j < repeats[i]; ++j) {
                blocks.push_back(patterns[i]);
            }
        }
        std::random_shuffle(blocks.begin(), blocks.end());
        for (const auto& block : blocks) {
            sequence += block;
        }
    }

    // Introduce mismatches
    std::srand(std::time(nullptr));
    for (int i = 0; i < mismatches; ++i) {
        int position = std::rand() % sequence.size();
        char original = sequence[position];
        char replacement;
        do {
            int randomNucleotide = std::rand() % 4;
            switch (randomNucleotide) {
                case 0: replacement = 'A'; break;
                case 1: replacement = 'T'; break;
                case 2: replacement = 'C'; break;
                case 3: replacement = 'G'; break;
            }
        } while (replacement == original);
        sequence[position] = replacement;
    }

    return sequence;
}

void generateFasta() {
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<std::string> patterns = askForPatterns();
    std::vector<int> repeats = askForRepeats(patterns.size());
    char mode = askForMode();
    int mismatches = askForMismatches();

    int totalLength = 0;
    for (size_t i = 0; i < patterns.size(); ++i) {
        totalLength += patterns[i].length() * repeats[i];
    }

    if (mismatches >= totalLength) {
        std::cout << "Mismatches are too many. Operation rejected." << std::endl;
        return;
    }

    std::string sequence = generateSequence(patterns, repeats, mode, mismatches);

    // Write to FASTA file
    std::ofstream fastaFile("/mnt/c/Users/xato_/Downloads/VGL/teloscope/testFiles/shortRandom.fasta");
    fastaFile << "> Telomeric patterns: ";
    for (const auto& pattern : patterns) {
        fastaFile << pattern << " ";
    }
    fastaFile << "| Number of repeats: ";
    for (const auto& repeat : repeats) {
        fastaFile << repeat << " ";
    }
    fastaFile << "| Mode: " << (mode == 'c' ? "concatenate" : "random") << "\n";
    fastaFile << sequence;  // Write the generated sequence
    fastaFile.close();

    auto finish = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

    std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;
}

int main() {
    generateFasta();
    return 0;
}