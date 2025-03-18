#include "main.h"
#include <input.h> // check
#include <iostream>


std::string version = "0.0.9";

// global
std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

int tabular_flag;
int verbose_flag;
int cmd_flag;

int maxThreads = 0;
std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;

Log lg;
std::vector<Log> logs;

UserInputTeloscope userInput; // init input object

int main(int argc, char **argv) {
    
    short int c; // optarg
    // short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    
    std::string cmd;

    bool isPipe = false; // to check if input is from pipe
    
    if (argc == 1) { // case: with no arguments
            
        printf("teloscope -f input.[fa/fa.gz] \nUse-h for additional help. \nUse -f to initiate the tool.");
        exit(0);
        
    }
    
    static struct option long_options[] = { // struct mapping long options
        {"input-sequence", required_argument, 0, 'f'},
        {"output", required_argument, 0, 'o'},
        {"patterns", required_argument, 0, 'p'},
        {"window", required_argument, 0, 'w'},
        {"step", required_argument, 0, 's'},
        {"canonical", required_argument, 0, 'c'},
        {"threads", required_argument, 0, 'j'},
        {"min-block-length", required_argument, 0, 'l'},
        {"max-block-distance", required_argument, 0, 'd'},
        {"terminal-limit", required_argument, 0, 't'},

        {"out-fasta", no_argument, 0, 'a'},
        {"out-win-repeats", no_argument, 0, 'r'},
        {"out-gc", no_argument, 0, 'g'},
        {"out-entropy", no_argument, 0, 'e'},
        {"out-matches", no_argument, 0, 'm'},
        {"out-its", no_argument, 0, 'i'},
        {"ultra-fast", no_argument, 0, 'u'},
        {"verbose", no_argument, &verbose_flag, 1},
        {"cmd", no_argument, &cmd_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    while (arguments) { // loop through argv
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "-:f:j:o:p:s:w:c:l:d:t:argemivhu", long_options, &option_index);
        
        if (c == -1) { // exit the loop if run out of options
            break;
            
        }
        switch (c) {
            case ':': // handle options without arguments
                switch (optopt) { // the command line option last matched
                    // case 'b':
                    //     break;
                        
                    default:
                        fprintf(stderr, "Error: Option -%c is missing a required argument\n", optopt);
                        return EXIT_FAILURE;
                }
                break;
            default: // handle positional arguments


            case 0: // case for long options without short options, none yet
                
                // if (strcmp(long_options[option_index].name,"line-length") == 0)
                //     splitLength = atoi(optarg);
                
                break;


            case 'f': // input sequence
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                    userInput.pipeType = 'f'; // pipe input is a sequence
                    
                } else { // input is a regular file
                    ifFileExists(optarg);
                    userInput.inSequence = optarg;
                }
                
                if (userInput.inSequence.empty()) {
                    fprintf(stderr, "Error: Input sequence file is required. Use -f or --input-sequence.\n");
                    exit(EXIT_FAILURE);
                }
                break;


            case 'o':
                {
                    userInput.outRoute = optarg;

                    if (userInput.outRoute.empty()) {
                        fprintf(stderr, "Error: Output route is required. Use -o or --output.\n");
                        exit(EXIT_FAILURE);
                    }

                    if (!std::filesystem::exists(userInput.outRoute)) {
                        std::filesystem::create_directories(userInput.outRoute); // create directory if it doesn't exist
                    }
                }
                break;


            case 'j': // max threads
                maxThreads = atoi(optarg);
                userInput.stats_flag = 1;
                break;


            case 'c': { // Handle canonical pattern
                if (!optarg || strlen(optarg) == 0) {
                    fprintf(stderr, "Warning: Empty canonical pattern provided, using default vertebrate TTAGGG.\n");
                } else {
                    std::string canonicalPattern = optarg;
                    unmaskSequence(canonicalPattern);

                    // Check for numerical characters
                    if (std::any_of(canonicalPattern.begin(), canonicalPattern.end(), ::isdigit)) {
                        fprintf(stderr, "Error: Canonical pattern '%s' contains numerical characters.\n", canonicalPattern.c_str());
                        exit(EXIT_FAILURE);
                    }

                    // Store canonical pattern and its reverse complement
                    userInput.canonicalSize = canonicalPattern.size();
                    userInput.canonicalFwd = canonicalPattern;
                    userInput.canonicalRev = revCom(canonicalPattern);
                    fprintf(stderr, "Setting canonical pattern: %s and its reverse complement: %s\n", userInput.canonicalFwd.c_str(), userInput.canonicalRev.c_str());
                }
                break;
            }

            
            case 'p': { // Handle search patterns
                if (!optarg || strlen(optarg) == 0) {
                    fprintf(stderr, "Warning: Empty pattern list provided, using default: TTAGGG, CCCTAA\n");
                
                } else {
                    userInput.patterns.clear(); // Clear default patterns first
                    
                    std::istringstream patternStream(optarg);
                    std::string pattern;
                    
                    while (std::getline(patternStream, pattern, ',')) {
                        if (pattern.empty()) continue;
            
                        if (std::any_of(pattern.begin(), pattern.end(), ::isdigit)) {
                            fprintf(stderr, "Error: Pattern '%s' contains numerical characters.\n", pattern.c_str());
                            exit(EXIT_FAILURE);
                        }
                        
                        unmaskSequence(pattern);
                        
                        // Generate all combinations for the pattern based on IUPAC codes
                        std::vector<std::string> combinations;
                        std::string current_pattern = pattern;
                        getCombinations(pattern, current_pattern, 0, combinations);
                        lg.verbose("Adding (" + std::to_string(combinations.size()) + ") telomeric patterns and their reverse complements");
            
                        // Add each combination and its reverse complement to userInput.patterns
                        for (const std::string &comb : combinations) {
                            userInput.patterns.emplace_back(comb);
                            userInput.patterns.emplace_back(revCom(comb));
                        }
                    }
            
                    if (!userInput.patterns.empty()) {
                        // Remove duplicates
                        std::sort(userInput.patterns.begin(), userInput.patterns.end());
                        auto last = std::unique(userInput.patterns.begin(), userInput.patterns.end());
                        userInput.patterns.erase(last, userInput.patterns.end());
                    }

                    // userInput.hammingDistances = getHammingDistances(userInput.patterns, userInput.canonicalFwd, userInput.canonicalRev); // To update
                    // lg.verbose("Hamming distances precomputed");
                }
                break;
            }


            case 'w': {
                try {
                    userInput.windowSize = std::stoi(optarg);
                    
                    if (userInput.windowSize <= 0) {
                        fprintf(stderr, "Error: Window size (-w or --window) must be > 0.\n");
                        exit(EXIT_FAILURE);
                    }
                } catch (const std::exception& e) {
                    fprintf(stderr, "Error: Invalid window size '%s'. Must be a number.\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            }


            case 's': {
                try {
                    userInput.step = std::stoi(optarg);
                    
                    if (userInput.step > userInput.windowSize) {
                        fprintf(stderr, "Error: Step size (%d) cannot be larger than window size (%d)!\n", userInput.step, userInput.windowSize);
                        exit(EXIT_FAILURE);
                        
                    } else if (userInput.step == userInput.windowSize) {
                        fprintf(stderr, "Warning: Equal step and window sizes will bin the sequence.\n");
                        fprintf(stderr, "Tip: Large sizes are recommended to avoid missing matches between windows.\n");
                        
                    } else {
                        fprintf(stderr, "Sliding with window size (%d) and step size (%d). \n", userInput.windowSize, userInput.step);
                    }
                    
                    if (userInput.step <= 0) {
                        fprintf(stderr, "Error: Step size (-s or --step) must be > 0.\n");
                        exit(EXIT_FAILURE);
                    }
                } catch (const std::exception& e) {
                    fprintf(stderr, "Error: Invalid step size '%s'. Must be a number.\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            }


            case 'l': {
                try {
                    userInput.minBlockLen = std::stoi(optarg);
                    
                    if (userInput.minBlockLen <= 0) {
                        fprintf(stderr, "Error: Min block length (-l or --min-block-length) must be > 0.\n");
                        exit(EXIT_FAILURE);
                    }
                } catch (const std::exception& e) {
                    fprintf(stderr, "Error: Invalid min block length '%s'. Must be a number.\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            }


            case 'd': {
                try {
                    userInput.maxBlockDist = std::stoi(optarg);
                    
                    if (userInput.maxBlockDist <= 0) {
                        fprintf(stderr, "Error: Max block distance (-d or --max-block-distance) must be > 0.\n");
                        exit(EXIT_FAILURE);
                    }
                } catch (const std::exception& e) {
                    fprintf(stderr, "Error: Invalid max block distance '%s'. Must be a number.\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            }


            case 't' : {
                try {
                    userInput.terminalLimit = std::stoi(optarg);
                    
                    if (userInput.terminalLimit <= 0) {
                        fprintf(stderr, "Error: Terminal limit (-t or --terminal-limit) must be > 0.\n");
                        exit(EXIT_FAILURE);
                    }
                } catch (const std::exception& e) {
                    fprintf(stderr, "Error: Invalid terminal limit '%s'. Must be a number.\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            }


            case 'a':
                userInput.outFasta = true;
                break;


            case 'r':
                userInput.outWinRepeats = true;
                break;


            case 'g':
                userInput.outGC = true;
                break;


            case 'e':
                userInput.outEntropy = true;
                break;


            case 'm':
                userInput.outMatches = true;
                break;


            case 'i':
                userInput.outITS = true;
                break;


            case 'u':
                userInput.ultraFastMode = true;
                userInput.outGC = false;
                userInput.outEntropy = false;
                userInput.outITS = false;
                userInput.outWinRepeats = false;
                userInput.outMatches = false;
                fprintf(stderr, "Warning: Ultra-fast mode enabled. Only scanning terminal regions.\n");
                fprintf(stderr, "Warning: Interstitial telomeres (ITS) scan and window-based metrics are disabled.\n");
                break;


            case 'v': // software version
                printf("/// Teloscope v%s\n", version.c_str());
                printf("\nDeveloped by:\nJack A. Medico amedico@rockefeller.edu\n");
                printf("\nDirected by:\nGiulio Formenti giulio.formenti@gmail.com\n");
                printf("\nhttps://www.vertebrategenomelab.org/home");
                exit(0);
                break;

            case 'h': // help
                printf("teloscope [commands]\n");
                printf("\nRequired Parameters:\n");
                printf("\t'-f'\t--input-sequence\tInitiate tool with fasta file.\n");
                printf("\t'-o'\t--output\tSet output route.\n");
                printf("\t'-c'\t--canonical\tSet canonical pattern. [Default: TTAGGG]\n");
                printf("\t'-p'\t--patterns\tSet patterns to explore, separate them by commas [Default: TTAGGG]\n");
                printf("\t'-w'\t--window\tSet sliding window size. [Default: 1000]\n");
                printf("\t'-s'\t--step\tSet sliding window step. [Default: 500]\n");
                printf("\t'-j'\t--threads\tSet maximum number of threads. [Default: max. available]\n");
                printf("\t'-l'\t--min-block-length\tSet minimum block length for merging. [Default: 500]\n");
                printf("\t'-d'\t--max-block-distance\tSet maximum block distance for merging. [Default: 50]\n");
                printf("\t'-t'\t--terminal-limit\tSet terminal limit for exploring telomere variant regions (TVRs). [Default: 50000]\n");

                printf("\nOptional Parameters:\n");
                printf("\t'-r'\t--out-win-repeats\tOutput canonical/noncanonical repeats and density by window. [Default: false]\n");
                printf("\t'-g'\t--out-gc\tOutput GC content for each window. [Default: false]\n");
                printf("\t'-e'\t--out-entropy\tOutput Shannon entropy for each window. [Default: false]\n");
                printf("\t'-m'\t--out-matches\tOutput all canonical and terminal non-canonical matches. [Default: false]\n");
                printf("\t'-i'\t--out-its\tOutput assembly interstitial telomere (ITSs) regions.[Default: false] \n");
                printf("\t'-u'\t--ultra-fast\tUltra-fast mode. Only scans terminal telomeres at contig ends. [Default: false]\n");

                printf("\t'-v'\t--version\tPrint current software version.\n");
                printf("\t'-h'\t--help\tPrint current software options.\n");
                printf("\t--verbose\tVerbose output.\n");
                printf("\t--cmd\tPrint command line.\n");
                exit(0);
        }
    }

    lg.verbose("Input variables assigned");
    
    if (cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }

    threadPool.init(maxThreads);

    Input in;
    in.load(userInput); // load user input
    lg.verbose("Loaded user input");
    
    InSequences inSequences; // initialize sequence collection object
    lg.verbose("Sequence object generated");
    in.read(inSequences); // read input content to inSequences container

    lg.verbose("Finished reading input files");
    if(verbose_flag) {std::cerr<<"\n";}; // giulio?

    threadPool.join();

    lg.verbose("Generated output");

    exit(EXIT_SUCCESS);
    
}