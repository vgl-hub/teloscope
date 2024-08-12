#include "main.h"
#include <input.h> // check


std::string version = "0.0.1";

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
std::string outRoute;

UserInputTeloscope userInput; // init input object

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    
    std::string cmd;

    bool isPipe = false; // to check if input is from pipe
    
    if (argc == 1) { // case: with no arguments
            
        printf("teloscope -f input.[fa/fa.gz] \n-h for additional help. Use -f to initiate the tool.\n");
        exit(0);
        
    }
    
    static struct option long_options[] = { // struct mapping long options
        {"input-sequence", required_argument, 0, 'f'},
        {"patterns", required_argument, 0, 'p'},
        {"window", required_argument, 0, 'w'},
        {"step", required_argument, 0, 's'},
        {"mode", required_argument, 0, 'm'},
        {"output", required_argument, 0, 'o'},
        {"threads", required_argument, 0, 'j'},
        
        {"verbose", no_argument, &verbose_flag, 1},
        {"cmd", no_argument, &cmd_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        // {"taxid", no_argument, 0, 't'},
        {0, 0, 0, 0}
    };
    
    while (arguments) { // loop through argv
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "-:f:j:m:o:p:s:w:vh", long_options, &option_index);

        // if (optind < argc && !isPipe) { // if pipe wasn't assigned already
            
        //     isPipe = isDash(argv[optind]) ? true : false; // check if the argument to the option is a '-' and set it as pipe input
            
        // }
        
        // if (optarg != nullptr && !isPipe) { // case where pipe input is given as positional argument (input sequence file)
        
        //     isPipe = isDash(optarg) ? true : false;
            
        // }
        
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

            case 0: // case for long options without short options, none yet
                
                // if (strcmp(long_options[option_index].name,"line-length") == 0)
                //     splitLength = atoi(optarg);
                
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
                std::set<std::string> providedModes;

                while (std::getline(modeStream, mode, ',')) {
                    if (mode.empty()) continue;

                    if (std::any_of(mode.begin(), mode.end(), ::isdigit)) {
                        std::cerr << "Error: Mode '" << mode << "' contains numerical characters.\n";
                        exit(EXIT_FAILURE);
                    }

                    std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);

                    if (mode == "all") {
                        providedModes = {"match", "entropy", "gc"};
                        break;
                    } else {
                        providedModes.insert(mode);
                    }
                }

                // Only disable modes that were not provided
                if (providedModes.find("match") == providedModes.end()) userInput.modeMatch = false;
                if (providedModes.find("entropy") == providedModes.end()) userInput.modeEntropy = false;
                if (providedModes.find("gc") == providedModes.end()) userInput.modeGC = false;

                break;
            }


            case 'o':
            {
                outRoute = optarg;
                
                if (outRoute.empty()) {
                    fprintf(stderr, "Error: Output route is required. Use --output or -o to specify it.\n"); // Jack: we have to define output as default
                    exit(EXIT_FAILURE);
                }
            }
                break;

            case 'p':
            {
                std::istringstream patternStream(optarg);
                std::string pattern;
                
                while (std::getline(patternStream, pattern, ',')) {
                    if (pattern.empty()) continue;

                    if (std::any_of(pattern.begin(), pattern.end(), ::isdigit)) {
                        std::cerr << "Error: Pattern '" << pattern << "' contains numerical characters.\n";
                        exit(EXIT_FAILURE);
                    }
                    
                    unmaskSequence(pattern);
                    
                    std::cout << "Adding pattern: " << pattern << " and its reverse complement" <<  "\n";
                    userInput.patterns.emplace_back(pattern);
                    userInput.patterns.emplace_back(revCom(pattern));
                }
                
                std::unique(userInput.patterns.begin(), userInput.patterns.end());
                if (userInput.patterns.empty()) {
                    userInput.patterns = {"TTAGGG", "CCCTAA"};
                    std::cout << "No patterns provided. Using canonical patterns: TTAGGG, CCCTAA" << "\n";
                }
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
                    fprintf(stderr, "Warning: Equal step and window sizes will bin the sequence.\n");
                    fprintf(stderr, "Tip: Large sizes are recommended to avoid missing matches between windows.\n");
                    
                } else {
                    fprintf(stderr, "Sliding windows with step size (%d) and window size (%d). \n", userInput.step, userInput.windowSize);
                    fprintf(stderr, "Tip: A step value close the window size results in faster runs.\n");
                }
                break;

            case 'v': // software version
                printf("/// Teloscope v%s\n", version.c_str());
                printf("Giulio Formenti giulio.formenti@gmail.com\n");
                printf("Jack A. Medico amedico@rockefeller.edu\n");
                exit(0);
                
            case 'h': // help
                printf("teloscope [commands]\n");
                printf("\nRequired Parameters:\n");
                printf("\t'-f'\t--input-sequence\tInitiate tool with fasta file.\n");
                printf("\t'-p'\t--patterns\tSet patterns to explore.\n");
                printf("\t'-w'\t--window\tSet sliding window size.\n");
                printf("\t'-s'\t--step\tSet sliding window step.\n");
                printf("\nOptional Parameters:\n");
                printf("\t'-v'\t--version\tPrint current software version.\n");
                printf("\t'-h'\t--help\tPrint current software options.\n");
                printf("\t--verbose\tverbose output.\n");
                exit(0);
        }
        
        if  (argc == 2 || // handle various cases in which the output should include summary stats
            (argc == 3 && pos_op == 2) ||
            (argc == 4 && pos_op == 3)) {
            
        }
        
    }

    lg.verbose("Input variables assigned");
    
    if (cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }

    threadPool.init(maxThreads); // initialize threadpool, giulio?

    Input in;
    in.load(userInput); // load user input
    lg.verbose("Loaded user input"); // jack: log not defined yet
    
    InSequences inSequences; // initialize sequence collection object
    lg.verbose("Sequence object generated");
    in.read(inSequences); // read input content to inSequences container

    lg.verbose("Finished reading input files");
    if(verbose_flag) {std::cerr<<"\n";}; // giulio?

    threadPool.join();

    lg.verbose("Generated output");

    exit(EXIT_SUCCESS);
    
}

