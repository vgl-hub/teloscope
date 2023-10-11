#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <filesystem>
#include <chrono>

struct TelomereInfo {
    std::string pattern;
    int repeats1, repeats2;
};

std::string generateRandomSequence(int length) { //, unsigned seed)
    std::string bases = "ATCG";
    std::string sequence;
    std::srand(std::time(nullptr)); //std::srand(seed);
    for (int i = 0; i < length; ++i) {
        sequence += bases[std::rand() % 4];
    }
    return sequence;
}

std::string introduceMismatches(const std::string& sequence, int mismatches) { //, unsigned seed)
    std::string mutatedSequence = sequence;
    std::srand(std::time(nullptr)); //std::srand(seed);
    for (int i = 0; i < mismatches; ++i) {
        int position = std::rand() % mutatedSequence.size();
        char original = mutatedSequence[position];
        char replacement;
        do {
            int randomNucleotide = std::rand() % 4;
            replacement = "ATCG"[randomNucleotide];
        } while (replacement == original);
        mutatedSequence[position] = replacement;
    }
    return mutatedSequence;
}

void readCSV(const std::string& fileName, int& chrSeqLength, char& mode1, char& mode2, double& mismatchRate, std::vector<TelomereInfo>& telomeres) {
    std::ifstream file(fileName); // functions.cpp gfalibs 366 > loop, store in a vector, change type later
    std::string line;
    std::getline(file, line);
    std::stringstream ss(line);
    std::string temp;
    std::getline(ss, temp, ','); chrSeqLength = std::stoi(temp);
    std::getline(ss, temp, ','); mode1 = temp[0];
    std::getline(ss, temp, ','); mode2 = temp[0];
    std::getline(ss, temp, ','); mismatchRate = std::stod(temp);
    
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        TelomereInfo info;
        std::getline(ss, info.pattern, ',');
        std::getline(ss, temp, ','); info.repeats1 = std::stoi(temp);
        std::getline(ss, temp, ','); info.repeats2 = std::stoi(temp);
        telomeres.push_back(info);
    }
}


std::string generateTelomereSequence(const TelomereInfo& info, char mode, int repeats) {
    std::string sequence;
    std::vector<std::string> blocks;

    if (mode == 'c') {
        for (int i = 0; i < repeats; ++i) {
            sequence += info.pattern;
        }
    }
    else if (mode == 'r') {
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
        
        for (int i = 0; i < repeats; ++i) {
            blocks.push_back(info.pattern);
        }
        
        std::random_shuffle(blocks.begin(), blocks.end());
        
        for (const auto& block : blocks) {
            sequence += block;
        }
    }
    return sequence;
}

void generateFasta(const std::string& fileName, const std::string& sequence, const std::string& headerInfo) {
    std::ofstream file(fileName); // get name, extracting from format 263 functions.cpp
    file << ">" << headerInfo << std::endl;
    file << sequence << std::endl;
}

int main() {
    //std::srand(static_cast<unsigned int>(std::time(nullptr)));
    // unsigned seed = std::time(nullptr);
    auto start = std::chrono::high_resolution_clock::now();


    int chrSeqLength;
    char mode1, mode2;
    double mismatchRate;
    std::vector<TelomereInfo> telomeres;
    
    readCSV("testfasta.csv", chrSeqLength, mode1, mode2, mismatchRate, telomeres);

    int tel1Length = 0, tel2Length = 0;
    std::string tel1, tel2;
    
    for (const auto& info : telomeres) {
        std::cout << "Debug: Processing telomere pattern = " << info.pattern << std::endl;
        tel1 += generateTelomereSequence(info, mode1, info.repeats1);
        tel2 += generateTelomereSequence(info, mode2, info.repeats2);
        tel1Length += info.pattern.length() * info.repeats1;
        tel2Length += info.pattern.length() * info.repeats2;
    }
    
    int randomSeqLength = chrSeqLength - tel1Length - tel2Length;
    
    if (randomSeqLength <= 0) {
        std::cerr << "Error: Non-telomeric sequence length must be positive.";
        return 1;
    }

    int mismatches1 = static_cast<int>(tel1Length * mismatchRate);
    int mismatches2 = static_cast<int>(tel2Length * mismatchRate);
    tel1 = introduceMismatches(tel1, mismatches1); //, seed);
    tel2 = introduceMismatches(tel2, mismatches2); //, seed);

    std::string randomSeq = generateRandomSequence(randomSeqLength); //, seed);
    std::string finalSeq = tel1 + randomSeq + tel2;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Debug: Number of telomere patterns = " << telomeres.size() << std::endl;
    std::cout << "Debug: tel1 = " << tel1 << std::endl;
    std::cout << "Debug: tel2 = " << tel2 << std::endl;
    std::cout << "Debug: tel1Length = " << tel1Length << std::endl;
    std::cout << "Debug: tel2Length = " << tel2Length << std::endl;
    std::cout << "Debug: randomSeqLength = " << randomSeqLength << std::endl;
    std::cout << "Debug: tel1 + mismatches = " << tel1 << std::endl;
    std::cout << "Debug: tel2 + mismatches = " << tel2 << std::endl;
    std::cout << "Debug: randomSeq = " << randomSeq << std::endl;
    std::cout << "Debug: finalSeq = " << finalSeq << std::endl; //add outputs in log?

    std::filesystem::path logPath = std::filesystem::current_path() / ".." / "testFiles" / "test.log";
    std::ofstream logFile(logPath.string(), std::ios_base::app);
    logFile << "Chromosome Length: " << chrSeqLength << ", Time Taken: " << elapsed.count() << "s\n";

    std::filesystem::path csvPath("testfasta.csv");
    std::string fastaFileName = csvPath.stem().string() + ".fasta";
    std::filesystem::path fastaPath = std::filesystem::current_path() / ".." / "testFiles" / fastaFileName;

    generateFasta(fastaPath.string(), finalSeq, "Header Information");

    return 0;
}
// add timestamps, large files