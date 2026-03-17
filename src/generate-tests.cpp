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
        {{"fasta", "fasta.gz", "fa", "fa.gz"}, {""}}
    //  {{set of test file extensions}, {list of command line args to run with}}
    };

    const std::map<std::set<std::string>, std::vector<std::string>> file_args = {
        {{"random1.fasta"}, {"-f testFiles/random1.fasta -w 3 -s 1"}},
        {{"random2.fasta"}, {"-f testFiles/random2.fasta -w 10 -s 5"}},
        {{"random3.fasta"}, {"-f testFiles/random3.fasta -w 10 -s 5"}},
        // Synthetic classification tests with flag combos
        {{"t2t.fa"}, {
            "-f testFiles/t2t.fa -w 500 -s 500",
            "-f testFiles/t2t.fa -w 500 -s 250",
            "-f testFiles/t2t.fa -w 100 -s 100",
            "-f testFiles/t2t.fa -l 1000",
            "-f testFiles/t2t.fa -l 100",
            "-f testFiles/t2t.fa -t 100",
            "-f testFiles/t2t.fa -t 3000",
            "-f testFiles/t2t.fa -k 1",
            "-f testFiles/t2t.fa -k 5",
            "-f testFiles/t2t.fa -k 10",
            "-f testFiles/t2t.fa -d 1",
            "-f testFiles/t2t.fa -d 5",
            "-f testFiles/t2t.fa -d 100",
            "-f testFiles/t2t.fa -n",
            "-f testFiles/t2t.fa -k 5 -d 5",
            "-f testFiles/t2t.fa -r -o testFiles/tmp",
            "-f testFiles/t2t.fa -w 200 -s 200 -r -o testFiles/tmp",
            "-f testFiles/t2t.fa -r -g -e -o testFiles/tmp",
            "-f testFiles/t2t.fa -i -o testFiles/tmp",
            "-f testFiles/t2t.fa -m -n -o testFiles/tmp",
            "-f testFiles/t2t.fa -y 0.99",
            "-f testFiles/t2t.fa -t 700 -l 100"
        }},
        {{"plant.fa"}, {
            "-f testFiles/plant.fa -c CCCTAAA",
            "-f testFiles/plant.fa -c CCCTAAA -p CCCTAAA",
            "-f testFiles/plant.fa -c CCCTAAA -p TTTAGGG,CCCTAAA"
        }},
        {{"edit_test.fa"}, {
            "-f testFiles/edit_test.fa -x 1",
            "-f testFiles/edit_test.fa -x 2",
            "-f testFiles/edit_test.fa -x 1 -p CCCTAA,TTAGGG"
        }},
        {{"density_edge.fa"}, {
            "-f testFiles/density_edge.fa -y 0.8",
            "-f testFiles/density_edge.fa -y 0.3",
            "-f testFiles/density_edge.fa -l 100 -y 0.3",
            "-f testFiles/density_edge.fa -l 1000 -y 0.3",
            "-f testFiles/density_edge.fa -d 1",
            "-f testFiles/density_edge.fa -l 500 -y 0.4 -k 10"
        }},
        {{"its.fa"}, {
            "-f testFiles/its.fa -t 100",
            "-f testFiles/its.fa -i -o testFiles/tmp",
            "-f testFiles/its.fa -m -i -o testFiles/tmp"
        }},
        {{"multi.fa"}, {
            "-f testFiles/multi.fa -n",
            "-f testFiles/multi.fa -l 100",
            "-f testFiles/multi.fa -r -o testFiles/tmp",
            "-f testFiles/multi.fa -k 5",
            "-f testFiles/multi.fa -d 1",
            "-f testFiles/multi.fa -m -o testFiles/tmp",
            "-f testFiles/multi.fa -i -o testFiles/tmp"
        }},
        {{"incomplete_p.fa"}, {
            "-f testFiles/incomplete_p.fa -k 1",
            "-f testFiles/incomplete_p.fa -n"
        }},
        {{"incomplete_q.fa"}, {
            "-f testFiles/incomplete_q.fa -m -o testFiles/tmp"
        }},
        {{"discordant.fa"}, {
            "-f testFiles/discordant.fa -t 3000",
            "-f testFiles/discordant.fa -n"
        }},
        {{"no_telo.fa"}, {
            "-f testFiles/no_telo.fa -m -o testFiles/tmp"
        }},
        {{"misassembly.fa"}, {
            "-f testFiles/misassembly.fa -m -o testFiles/tmp"
        }},
        {{"short_contig.fa"}, {
            "-f testFiles/short_contig.fa -w 50 -s 50"
        }},
        {{"balanced.fa"}, {
            "-f testFiles/balanced.fa -i -o testFiles/tmp"
        }},
        {{"gapped_t2t.fa"}, {
            "-f testFiles/gapped_t2t.fa -r -o testFiles/tmp",
            "-f testFiles/gapped_t2t.fa -l 100",
            "-f testFiles/gapped_t2t.fa -n",
            "-f testFiles/gapped_t2t.fa -r -g -e -o testFiles/tmp"
        }},
        // Directional scanning: mirrored/inverted p/q edge cases
        {{"mirror_inverted_short.fa"}, {
            "-f testFiles/mirror_inverted_short.fa -n",
            "-f testFiles/mirror_inverted_short.fa -r -o testFiles/tmp"
        }},
        {{"mirror_rev_start.fa"}, {
            "-f testFiles/mirror_rev_start.fa -n",
            "-f testFiles/mirror_rev_start.fa -r -o testFiles/tmp"
        }},
        {{"mirror_inverted_long.fa"}, {
            "-f testFiles/mirror_inverted_long.fa -t 1000 -n",
            "-f testFiles/mirror_inverted_long.fa -t 1000 -r -o testFiles/tmp"
        }},
        {{"mirror_fwd_end_long.fa"}, {
            "-f testFiles/mirror_fwd_end_long.fa -t 1000 -n"
        }},
        {{"mirror_rev_start_long.fa"}, {
            "-f testFiles/mirror_rev_start_long.fa -t 1000 -n"
        }},
        {{"mirror_both_start.fa"}, {
            "-f testFiles/mirror_both_start.fa -n",
            "-f testFiles/mirror_both_start.fa -r -o testFiles/tmp"
        }},
        {{"mirror_both_end.fa"}, {
            "-f testFiles/mirror_both_end.fa -n",
            "-f testFiles/mirror_both_end.fa -r -o testFiles/tmp"
        }},
        {{"mirror_extend.fa"}, {
            "-f testFiles/mirror_extend.fa -t 300 -n",
            "-f testFiles/mirror_extend.fa -t 300 -r -o testFiles/tmp"
        }},
        {{"mirror_merge.fa"}, {
            "-f testFiles/mirror_merge.fa -n",
            "-f testFiles/mirror_merge.fa -r -o testFiles/tmp"
        }},
        {{"mirror_no_merge.fa"}, {
            "-f testFiles/mirror_no_merge.fa -n"
        }}
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
