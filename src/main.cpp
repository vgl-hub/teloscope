#include <input.h>
#include <main.h>

std::string version = "0.1.2";

std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
short int tabular_flag; // jack: why don't we group this by struct (e.g. input, output, operations)?
int verbose_flag;
int seqReport_flag;
int outSequence_flag;
int nstarReport_flag;
int outSize_flag;
int outCoord_flag;
int outFile_flag;
int outBubbles_flag;
int stats_flag;
int rmGaps_flag;
int discoverPaths_flag;
int extractContigs_flag;
int hc_flag;
int hc_cutoff;
int terminalOvlLen = 0;
int maxThreads = 0;
int cmd_flag;

std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;
Log lg;

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    
    std::string cmd;

    bool isPipe = false; // to check if input is from pipe

    UserInputTeloscope userInput; //UserInputTeloscope
    
    if (argc == 1) { // mytool with no arguments
            
        printf("teloscope [command]\n-h for additional help. Use -f to initiate the tool.\n");
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
        // {"taxid", no_argument, 0, 't'},

        {0, 0, 0, 0}
    };
    
    while (arguments) { // loop through argv
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "-:f:pw:s:vh", long_options, &option_index);
        
        if (c == -1) { // exit the loop if run out of options
            break;
            
        }

        switch (c) {
            case ':': // handle options without arguments
                switch (optopt) { // the command line option last matched
                    case 'b':
                        break;
                        
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
                
            case 'p': {
                // User defined-patterns
                std::istringstream patternStream(optarg);
                std::string pattern;
                while (std::getline(patternStream, pattern, ',')) {
                    userInput.patterns.push_back(pattern);
                }
                break;
            }

            case 'w':
                userInput.windowSize = std::stoi(optarg);
                break;

            case 's':
                userInput.step = std::stoi(optarg);
                break;
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

