#include <input.h> // check
#include <main.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>

std::string version = "0.1.3";
std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

short int tabular_flag;
int verbose_flag;
int cmd_flag;

int maxThreads = 0;
std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;
Log lg;

UserInputTeloscope userInput; // init input object

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    
    std::string cmd;

    bool isPipe = false; // to check if input is from pipe
    
    if (argc == 1) { // case: with no arguments
            
        printf("teloscope -f input.[fasta] \n-h for additional help. Use -f to initiate the tool.\n");
        exit(0);
        
    }
    
    static struct option long_options[] = { // struct mapping long options
        {"input-sequence", required_argument, 0, 'f'},
        {"verbose", no_argument, &verbose_flag, 1},
        {"cmd", no_argument, &cmd_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'}, // giulio: expand the arguments, tax ID
        {"patterns", required_argument, 0, 'p'},
        {"window", required_argument, 0, 'w'},
        {"step", required_argument, 0, 's'},
        {"mode", required_argument, 0, 'm'},
        {"output", required_argument, 0, 'o'},
        {"threads", required_argument, 0, 'j'},
        // {"taxid", no_argument, 0, 't'},
        {0, 0, 0, 0}
    };
    
    while (arguments) { // loop through argv
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "-:f:j:m:o:p:s:w:vh", long_options, &option_index);
        
        if (c == -1) { // exit the loop if run out of options
            break;
            
        }
        switch (c) {
            case ':': // handle options without arguments
                switch (optopt) { // the command line option last matched
                    // case 'b':
                    //     break;
                        
                    default:
                        fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                        return EXIT_FAILURE;
                }
                break;
            default: // handle positional arguments
                                
            case 0: // case for long options without short options
                
//                if (strcmp(long_options[option_index].name,"line-length") == 0)
//                  splitLength = atoi(optarg);
                
                break;

            case 'f': // input sequence
                
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'f'; // pipe input is a sequence
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    userInput.inSequence = optarg;
                    
                }
                    
                break;
            
            case 'j': // max threads
                maxThreads = atoi(optarg);
                userInput.stats_flag = 1;
                break;

            case 'm': {
                std::istringstream modeStream(optarg);
                std::string mode;
                bool allModeSelected = false;
                while (std::getline(modeStream, mode, ',')) {
                    if (mode == "all") {
                        allModeSelected = true;
                        break;
                    }
                    userInput.mode.push_back(mode);
                }
                if (allModeSelected || userInput.mode.empty()) {
                    userInput.mode = {"match", "count", "density", "entropy", "gc"};
                }
                break;
            }

            case 'o':
            {
                std::string outRoute;
                outRoute = optarg;
                
                if (outRoute.empty()) {
                    fprintf(stderr, "Error: Output route is required. Use --output or -o to specify it.\n");
                    exit(EXIT_FAILURE);
                }
            }
                break;
                
            case 'p':
            {
                std::istringstream patternStream(optarg);
                std::string pattern;
                while (std::getline(patternStream, pattern, ',')) {
                    std::cout << "Adding pattern: " << pattern << std::endl;
                    userInput.patterns.emplace_back(pattern);
                }
                // if (userInput.patterns.empty()) {
                //     userInput.patterns = {"TTAGGG", "CCCTAA"};
                //     std::cout << "No patterns selected: Using canonical TTAAGGG and CCCTAA" << std::endl;
                // }
            }
                break;

            case 'w':
                userInput.windowSize = std::stoi(optarg);
                break;

            case 's':
                userInput.step = std::stoi(optarg);
                printf("/// Teloscope v%s\n", version.c_str());

                if (userInput.step > userInput.windowSize) {
                    fprintf(stderr, "Error: Step size (%d) larger than window size (%d) is not allowed.\n", userInput.step, userInput.windowSize);
                    exit(EXIT_FAILURE);

                } else if (userInput.step == userInput.windowSize) {
                    fprintf(stderr, "Warning: Equal step and window sizes will bin the sequence. Large window/step sizes are recommended to avoid missing matches.\n");
                    
                } else {
                    fprintf(stderr, "Sliding windows with step size (%d) and window size (%d). \n", userInput.step, userInput.windowSize);
                    fprintf(stderr, "Note: A step value close the window size results in fast runs.\n");
                }
                break;

            case 'v': // software version
                printf("/// Teloscope v%s\n", version.c_str());
                printf("Giulio Formenti giulio.formenti@gmail.com\n");
                printf("Jack A. Medico amedico@rockefeller.edu\n");
                exit(0);
                
            case 'h': // help
                printf("teloscope [command]\n");
                printf("\nRequired Parameters:\n");
                printf("--input-sequence    '-f'    Initiate tool with fasta file.\n");
                printf("--patterns    '-p'    Set patterns to explore.\n");
                printf("--window    '-w'    Set sliding window size.\n");
                printf("--step    '-s'    Set sliding window step.\n");
                printf("\nOptional Parameters:\n");
                printf("--version   -v  Print current software version.\n");
                printf("--help   -h  Print current software options.\n");
                exit(0);
        }
        
        if  (argc == 2 || // handle various cases in which the output should include summary stats
            (argc == 3 && pos_op == 2) ||
            (argc == 4 && pos_op == 3)) {
            
        }
        
    }
    
    if (cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }

    threadPool.init(maxThreads); // initialize threadpool

    Input in;
    in.load(userInput); // load user input
    lg.verbose("Loaded user input");
    
    InSequences inSequences; // initialize sequence collection object
    lg.verbose("Sequence object generated");
    in.read(inSequences); // read input content to inSequences container

    threadPool.join();

    exit(EXIT_SUCCESS);
    
}

