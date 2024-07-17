#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <random>
#include <thread>
#include <cmath>
#include <queue>
#include "../gfalibs/include/log.h"

struct TelomereInfo {
    std::string pattern;
    int pRepeats, qRepeats;
};

Log lg;

std::string getRandomSequence(int length);
void readCSV(const std::string& fileName, std::vector<TelomereInfo>& telomeres);
std::string generateAndMutateTelomere(const std::vector<TelomereInfo>& telomeres, char mode, bool isPTelomere, double substitutionRate);


int main(int argc, char **argv) {
    auto start = std::chrono::high_resolution_clock::now();
    
    ThreadPool<std::function<void()>> pool;  // type??
    pool.init(std::thread::hardware_concurrency());  // Initialize the thread pool

    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "print help")
        ("input,i", po::value<std::string>()->required(), "Input CSV file")
        ("length,l", po::value<int>()->required(), "Chromosome length")
        ("mode,m", po::value<std::string>()->default_value("concatenate"), "Mode ('concatenate' or 'random')")
        ("substitution,s", po::value<double>()->default_value(0.01), "Substitution rate")
        ("output,o", po::value<std::string>(), "Output FASTA file");
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    std::string inputFileName = vm["input"].as<std::string>();
    int chromosomeLength = vm["length"].as<int>();
    std::string mode = vm["mode"].as<std::string>();
    double substitutionRate = vm["substitution"].as<double>();
    std::string outputFileName = vm.count("output") ? vm["output"].as<std::string>() : inputFileName + ".fasta";

    std::vector<TelomereInfo> telomeres;
    readCSV(inputFileName, telomeres);

    std::string pTelomere, qTelomere;

    int nonTelomereLength = chromosomeLength - pTelomere.length() - qTelomere.length();
    std::string nonTelomere = getRandomSequence(nonTelomereLength);
    nonTelomere = addSubstitutions(nonTelomere, substitutionRate);

    std::string finalSeq = pTelomere + nonTelomere + qTelomere;

    std::ofstream fastaFile("../testFiles/" + outputFileName);
    fastaFile << "> Length: " << chromosomeLength << ", mode: " << mode << ", substitution rate: " << substitutionRate << "\n";
    for (size_t i = 0; i < finalSeq.length(); i += 80) {
        fastaFile << finalSeq.substr(i, 80) << "\n";
    }
    fastaFile.close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::ofstream logFile("../testFiles/test.log", std::ios_base::app);
    logFile << elapsed_seconds.count() << "\t" << chromosomeLength << "\t"
            << (pTelomere.length() + qTelomere.length()) << "\t"
            << ((double)chromosomeLength / elapsed_seconds.count()) << "bp/ms"
            << std::endl;
    logFile.close();

    return 0;
}

std::string getRandomSequence(int length) {
    static std::mt19937 gen(std::random_device{}());
    static std::uniform_int_distribution<> dist(0, 3);
    const std::string bases = "ATCG";
    std::string sequence;
    for (int i = 0; i < length; ++i) {
        sequence += bases[dist(gen)];
    }
    return sequence;
}


void readCSV(const std::string& fileName, std::vector<TelomereInfo>& telomeres) {
    std::ifstream file(fileName);
    if (!file) {
        throw std::runtime_error("Failed to open file: " + fileName);
    }
    std::string line;
    std::getline(file, line);  // Skip the header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        TelomereInfo info;
        std::string temp;
        std::getline(ss, info.pattern, ',');
        std::getline(ss, temp, ','); info.pRepeats = std::stoi(temp);
        std::getline(ss, temp, ','); info.qRepeats = std::stoi(temp);
        telomeres.push_back(info);
    }
}

std::string generateAndMutateTelomere(const std::vector<TelomereInfo>& telomeres, char mode, bool isPTelomere, double substitutionRate) {
    std::string sequence;
    std::vector<std::string> blocks;

    // Generate telomere
    for (const auto& info : telomeres) {
        int repeats = isPTelomere ? info.pRepeats : info.qRepeats;
        for (int i = 0; i < repeats; ++i) {
            blocks.push_back(info.pattern);
        }
    }

    if (mode == 'r') {
        std::shuffle(blocks.begin(), blocks.end(), std::mt19937{std::random_device{}()});
    }

    for (const auto& block : blocks) {
        sequence += block;
    }

    // Add substitutions
    static std::mt19937 gen(std::random_device{}());
    int substitutions = std::round(sequence.length() * substitutionRate);
    static std::uniform_int_distribution<> posDist(0, sequence.length() - 1);
    static std::uniform_int_distribution<> baseDist(0, 3);
    const std::string bases = "ATCG";

    for (int i = 0; i < substitutions; ++i) {
        int position = posDist(gen);
        char replacement;
        do {
            replacement = bases[baseDist(gen)];
        } while (replacement == sequence[position]);
        sequence[position] = replacement;
    }

    return sequence;
}
