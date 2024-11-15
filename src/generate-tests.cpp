#include <stdlib.h>
#include <map>
#include <cstdio>

#include "validate.h"

int main(int, char **argv) {
    std::cout << "WARNING: only run this program if teloscope is in a working state" << std::endl;
    std::cout << "WARNING: previous validate files will be deleted" << std::endl;
    std::cout << "continue? (Y/N) ";
    std::string input;
    std::cin >> input;
    if(input != "Y" && input != "y") {
        std::cout << "validate generation cancelled" << std::endl;
        std::exit(0);
    }

    std::cout << "deleting old validate files..." << std::endl;
    for(auto &file : list_dir("validateFiles")) {
        if(getFileExt(file) != "tst") continue; // dont delete README
        file = "validateFiles/"+file;
        if(remove(file.c_str()) != 0) {
            std::cerr << "error deleting <" << file << ">" << std::endl;
            return -1;
        }
    }

    std::cout << "generating new validate files..." << std::endl;

    std::string exePath = getExePath(argv[0]);

    const std::map<std::set<std::string>, std::vector<std::string>> ext_args = {
        {{"fasta", "fasta.gz"}, {""}}
    //  {{set of test file extensions}, {list of command line args to run with}}
    };

    const std::map<std::set<std::string>, std::vector<std::string>> file_args = {
        {{"random1.fasta"}, {"-f testFiles/random1.fasta -w 3 -s 1"}},
        {{"random2.fasta"}, {"-f testFiles/random2.fasta -w 10 -s 5"}},
        {{"random3.fasta"}, {"-f testFiles/random3.fasta -w 10 -s 5"}}
    //  {{set of test file paths}, {list of command line args to run with}}
    };

    const std::set<std::string> exclude {"agp", "sak"};

    for(const std::string &file : list_dir("testFiles")) {
        std::string ext = getFileExt(file);
        if(exclude.count(ext)) continue;
        for(auto pair : ext_args) {
            if(!pair.first.count(ext)) continue;
            for(auto args : pair.second) {
                genTest(exePath, file, args);
            }
        }
    }

    std::fstream fstream;
    for(const auto &pair : file_args) {
        for(const std::string &file : pair.first) {
            fstream.open("testFiles/"+file);
            if(!fstream) continue;
            fstream.close();
            for(const std::string &args : pair.second) {
                genTest(exePath, file, args);
            }
        }
    }

    std::exit(EXIT_SUCCESS);
}




// // JACK: New functions
// void cleanOldTests() {
//     std::vector<std::string> oldTests = {"validateFiles/random1.fasta.1.tst", "validateFiles/random2.fasta.1.tst"};
//     for (const auto& testFile : oldTests) {
//         if (remove(testFile.c_str()) != 0) {
//             std::cerr << "Error deleting " << testFile << std::endl;
//         }
//     }
// }
// void generateTest(const std::string& exePath, const std::string& inputFile, const std::string& args) {
//     std::string tstFile = "validateFiles/" + inputFile + ".tst";
//     std::ofstream out(tstFile);
//     if (!out) {
//         std::cerr << "Failed to create test file: " << tstFile << std::endl;
//         return;
//     }
//     out << exePath << " testFiles/" << inputFile << " " << args << " > embedded\n";
//     out.close();
//     std::cout << "Generated test file: " << tstFile << std::endl;
// }
// int main(int argc, char** argv) {
//     if (argc < 2) {
//         std::cerr << "Usage: " << argv[0] << " <path_to_teloscope_executable>" << std::endl;
//         return EXIT_FAILURE;
//     }
//     std::string exePath = argv[1];
//     std::cout << "WARNING: Previous validate files will be deleted. Continue? (Y/N) ";
//     std::string input;
//     std::cin >> input;
//     if (input != "Y" && input != "y") {
//         std::cout << "Validate generation cancelled." << std::endl;
//         return EXIT_FAILURE;
//     }
//     cleanOldTests();
//     std::cout << "Generating new validate files..." << std::endl;
//     // Generate tests based on specific command line arguments
//     generateTest(exePath, "random1.fasta", "-w 3 -s 1");
//     generateTest(exePath, "random2.fasta", "-p TTAGGG,CCCTAA -w 10 -s 5");
//     return EXIT_SUCCESS;
// }
