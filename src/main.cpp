#include "main.h"
#include <input.h>
#include <iostream>

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

std::string version = "0.1.4";

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

static std::string findReportScript(const char* argv0) {
    std::filesystem::path exeDir;

#if defined(__APPLE__)
    {
        char buf[4096];
        uint32_t bufsize = sizeof(buf);
        if (_NSGetExecutablePath(buf, &bufsize) == 0) {
            std::error_code ec;
            auto p = std::filesystem::canonical(buf, ec);
            if (!ec) exeDir = p.parent_path();
        }
    }
#elif !defined(_WIN32)
    {
        std::error_code ec;
        auto p = std::filesystem::canonical("/proc/self/exe", ec);
        if (!ec) exeDir = p.parent_path();
    }
#endif

    if (exeDir.empty() && argv0) {
        std::error_code ec;
        auto p = std::filesystem::canonical(argv0, ec);
        if (!ec) exeDir = p.parent_path();
    }

    if (exeDir.empty()) return "";

    for (const char* rel : {"../../scripts/teloscope_report.py",
                            "../scripts/teloscope_report.py",
                            "scripts/teloscope_report.py",
                            "teloscope_report.py"}) {
        std::error_code ec;
        auto resolved = std::filesystem::canonical(exeDir / rel, ec);
        if (!ec) return resolved.string();
    }

    return "";
}

int main(int argc, char **argv) {
    
    short int c; // optarg
    // short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    
    std::string cmd;

#ifdef _WIN32
    bool isPipe = !_isatty(_fileno(stdin));
#else
    bool isPipe = !isatty(STDIN_FILENO);
#endif
    bool hasInputPatterns = false;
    
    if (argc == 1 && !isPipe) { // case: with no arguments and no pipe

        printf("teloscope input.[fa|fa.gz|gfa] [options]\nteloscope -f input.[fa|fa.gz|gfa] -o /output/path/\nUse -h for additional help.\n");
        exit(0);

    }

    auto setInputFile = [&](const char* path) {
        ifFileExists(path);
        userInput.inSequence = path;
        const std::filesystem::path real = std::filesystem::canonical(userInput.inSequence);
        userInput.inSequence       = real.string();
        userInput.inSequencePrefix = real.parent_path().string();
        userInput.inSequenceName   = real.filename().string();
        userInput.outRoute         = userInput.inSequencePrefix;
    };
    
    static struct option long_options[] = { // struct mapping long options
        {"input-sequence", required_argument, 0, 'f'},
        {"output", required_argument, 0, 'o'},
        {"patterns", required_argument, 0, 'p'},
        {"window", required_argument, 0, 'w'},
        {"step", required_argument, 0, 's'},
        {"canonical", required_argument, 0, 'c'},
        {"threads", required_argument, 0, 'j'},
        {"terminal-limit", required_argument, 0, 't'},
        {"max-match-distance", required_argument, 0, 'k'},
        {"max-block-distance", required_argument, 0, 'd'},
        {"min-block-length", required_argument, 0, 'l'},
        {"min-block-density", required_argument, 0, 'y'},
        {"edit-distance", required_argument, 0, 'x'},

        {"out-fasta", no_argument, 0, 'a'},
        {"out-win-repeats", no_argument, 0, 'r'},
        {"out-gc", no_argument, 0, 'g'},
        {"out-entropy", no_argument, 0, 'e'},
        {"out-matches", no_argument, 0, 'm'},
        {"out-its", no_argument, 0, 'i'},
        {"ultra-fast", no_argument, 0, 'u'},
        {"manual-curation", no_argument, 0, 'n'},
        {"plot-report", no_argument, 0, 0},
        {"verbose", no_argument, &verbose_flag, 1},
        {"cmd", no_argument, &cmd_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    while (arguments) { // loop through argv
        
        int option_index = 0;

        c = getopt_long(argc, argv, "-:f:j:o:p:s:w:c:t:k:d:l:y:x:argemivhun", long_options, &option_index);

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


            case 0: // long options without short options
                if (strcmp(long_options[option_index].name, "plot-report") == 0)
                    userInput.outPlotReport = true;
                break;


            case 1: // positional argument (non-option)
                if (userInput.inSequence.empty()) {
                    setInputFile(optarg);
                } else {
                    fprintf(stderr, "Warning: Ignoring extra positional argument '%s'\n", optarg);
                }
                break;


            case 'f': // input sequence
                setInputFile(optarg);

                if (userInput.inSequence.empty()) {
                    fprintf(stderr, "Error: Input sequence file is required. Use -f or --input-sequence.\n");
                    exit(EXIT_FAILURE);
                }
                break;


            case 'o': // output route
                {
                    userInput.outRoute = optarg;

                    if (userInput.outRoute.empty()) {
                        userInput.outRoute = userInput.inSequencePrefix;
                        exit(EXIT_FAILURE);
                    }

                    try {
                        if (!std::filesystem::exists(userInput.outRoute)) {
                            std::filesystem::create_directories(userInput.outRoute);
                        }
                    } catch (const std::filesystem::filesystem_error& e) {
                        fprintf(stderr, "Error: Cannot create output directory '%s': %s\n",
                                userInput.outRoute.c_str(), e.code().message().c_str());
                        exit(EXIT_FAILURE);
                    }
                }
                break;


            case 'j': // max threads
                maxThreads = atoi(optarg);
                userInput.stats_flag = 1;
                break;


            case 'c': { // canonical pattern
                if (!optarg || strlen(optarg) == 0) {
                    fprintf(stderr, "Warning: Empty canonical pattern provided, using default vertebrate TTAGGG.\n");
                } else {
                    std::string canonicalPattern = optarg;
                    unmaskSequence(canonicalPattern);

                    if (std::any_of(canonicalPattern.begin(), canonicalPattern.end(), ::isdigit)) {
                        fprintf(stderr, "Error: Canonical pattern '%s' contains numerical characters.\n", canonicalPattern.c_str());
                        exit(EXIT_FAILURE);
                    }

                    // lex-smaller = Fwd
                    userInput.canonicalSize = canonicalPattern.size();
                    std::string revComp = revCom(canonicalPattern);
                    if (canonicalPattern <= revComp) {
                        userInput.canonicalFwd = canonicalPattern;
                        userInput.canonicalRev = revComp;
                    } else {
                        userInput.canonicalFwd = revComp;
                        userInput.canonicalRev = canonicalPattern;
                    }
                    fprintf(stderr, "Setting canonical pattern: %s and its reverse complement: %s\n", userInput.canonicalFwd.c_str(), userInput.canonicalRev.c_str());
                }
                break;
            }

            
            case 'p': { // search patterns
                hasInputPatterns = true;
                if (!optarg || strlen(optarg) == 0) {
                    fprintf(stderr, "Warning: Empty pattern list provided, using default: TTAGGG, CCCTAA\n");
                } else {
                    userInput.rawPatterns.clear();
                    std::istringstream patternStream(optarg);
                    std::string pattern;

                    while (std::getline(patternStream, pattern, ',')) {
                        if (pattern.empty()) continue;

                        if (std::any_of(pattern.begin(), pattern.end(), ::isdigit)) {
                            fprintf(stderr, "Error: Pattern '%s' contains numerical characters.\n", pattern.c_str());
                            exit(EXIT_FAILURE);
                        }

                        unmaskSequence(pattern);
                        userInput.rawPatterns.emplace_back(pattern);
                    }
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


            case 'k': { // max match distance
                try {
                    int v = std::stoi(optarg);
                    if (v <= 0) {
                        fprintf(stderr, "Error: Max match distance (-k/--max-match-distance) must be > 0.\n");
                        exit(EXIT_FAILURE);
                    }
                    userInput.maxMatchDist = static_cast<unsigned short>(v);
                } catch (...) {
                    fprintf(stderr, "Error: Invalid max match distance '%s'. Must be a number.\n", optarg);
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


            case 'y': { // min block density
                try {
                    float v = std::stof(optarg);
                    if (v < 0.0f || v > 1.0f) {
                        fprintf(stderr, "Error: Min block density (-y/--min-block-density) must be in the range [0,1].\n");
                        exit(EXIT_FAILURE);
                    }
                    userInput.minBlockDensity = v;
                } catch (...) {
                    fprintf(stderr, "Error: Invalid min block density '%s'. Must be a number [0,1].\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            }


            case 'x': { // edit distance
                try {
                    int v = std::stoi(optarg);
                    if (v < 0 || v > 2) {
                        fprintf(stderr, "Error: Edit distance (-x/--edit-distance) must be in the range [0,2].\n");
                        exit(EXIT_FAILURE);
                    }
                    userInput.editDistance = static_cast<uint8_t>(v);
                } catch (...) {
                    fprintf(stderr, "Error: Invalid edit distance '%s'. Must be a number [0,2].\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            }


            case 'a':
                userInput.outFasta = true;
                userInput.ultraFastMode = false;
                break;


            case 'r':
                userInput.outWinRepeats = true;
                userInput.ultraFastMode = false;
                break;


            case 'g':
                userInput.outGC = true;
                userInput.ultraFastMode = false;
                break;


            case 'e':
                userInput.outEntropy = true;
                userInput.ultraFastMode = false;
                break;


            case 'm':
                userInput.outMatches = true;
                userInput.ultraFastMode = false;
                break;


            case 'i':
                userInput.outITS = true;
                userInput.ultraFastMode = false;
                break;


            case 'u': {
                if (userInput.outWinRepeats || userInput.outGC ||
                    userInput.outEntropy   || userInput.outITS ||
                    userInput.outMatches) {
                    // conflicts with genome-wide flags, ignore -u
                    userInput.ultraFastMode = false;
                    fprintf(stderr, "Ignoring -u: -r/-g/-e/-i/-m request genome-wide scanning.\n");
                } else {
                    // terminal-only mode
                    userInput.ultraFastMode = true;
                    fprintf(stderr, "Fast mode: Only scanning terminal regions.\n");
                }
                break;
            }


            case 'n': // manual curation mode
                userInput.manualCuration = true;
                break;


            case 'v': // software version
                printf("/// Teloscope v%s\n", version.c_str());
                printf("\nDeveloped by:\nJack A. Medico amedico@rockefeller.edu\n");
                printf("\nDirected by:\nGiulio Formenti giulio.formenti@gmail.com\n");
                printf("\nhttps://www.vertebrategenomelab.org/home");
                exit(0);
                break;

            case 'h': // help
                printf("teloscope input.[fa|fa.gz|gfa] [options]\n");
                printf("teloscope -f input.[fa|fa.gz|gfa] [options]\n");
                printf("\nRequired Parameters:\n");
                printf("\t'-f'\t--input-sequence\tInput fasta file (or pass as first positional argument).\n");
                printf("\t'-o'\t--output\tSet output route. [Default: Input path]\n");
                printf("\t'-c'\t--canonical\tSet canonical pattern. [Default: TTAGGG]\n");
                printf("\t'-p'\t--patterns\tSet patterns to explore, separate them by commas [Default: TTAGGG]\n");
                printf("\t'-j'\t--threads\tSet maximum number of threads. [Default: max. available]\n");
                printf("\t'-t'\t--terminal-limit\tSet terminal limit for exploring telomere variant regions (TVRs). [Default: 50000]\n");
                printf("\t'-k'\t--max-match-distance\tSet maximum distance for merging matches. [Default: 50]\n");
                printf("\t'-d'\t--max-block-distance\tSet maximum block distance for extension. [Default: 200]\n");
                printf("\t'-l'\t--min-block-length\tSet minimum block length. [Default: 500]\n");
                printf("\t'-y'\t--min-block-density\tSet minimum block density. [Default: 0.5]\n");
                printf("\t'-x'\t--edit-distance\tSet edit distance for pattern matching (0-2). [Default: 1]\n");

                printf("\nOptional Parameters:\n");
                printf("\t'-w'\t--window\tSet sliding window size. [Default: 1000]\n");
                printf("\t'-s'\t--step\tSet sliding window step. [Default: 1000 (non-overlapping)]\n");
                printf("\t'-r'\t--out-win-repeats\tOutput per-window repeat density, canonical ratio, and strand ratio. [Default: false]\n");
                printf("\t'-g'\t--out-gc\tOutput GC content for each window. [Default: false]\n");
                printf("\t'-e'\t--out-entropy\tOutput Shannon entropy for each window. [Default: false]\n");
                printf("\t'-m'\t--out-matches\tOutput all canonical and terminal non-canonical matches. [Default: false]\n");
                printf("\t'-i'\t--out-its\tOutput assembly interstitial telomere (ITSs) regions.[Default: false] \n");
                printf("\t'-u'\t--ultra-fast\tUltra-fast mode. Only scans terminal telomeres at contig ends. [Default: true]\n");
                printf("\t'-n'\t--manual-curation\tRetain all terminal telomeres (contig + scaffold) in BED output. [Default: scaffold only]\n");
                printf("\t\t--plot-report\tGenerate a PDF plot report after analysis (requires Python 3 + matplotlib). [Default: false]\n");

                printf("\t'-v'\t--version\tPrint current software version.\n");
                printf("\t'-h'\t--help\tPrint current software options.\n");
                printf("\t--verbose\tVerbose output.\n");
                printf("\t--cmd\tPrint command line.\n");
                exit(0);
        }
    }

    // pipe-only invocation
    if (userInput.inSequence.empty() && isPipe) {
        userInput.pipeType = 'f';
        userInput.inSequenceName = "stdin";
        if (userInput.outRoute.empty())
            userInput.outRoute = ".";
    }

    // gzipped stdin not supported
    if (isPipe && userInput.pipeType == 'f') {
        int b = std::cin.peek();
        if (b == 0x1f) {
            fprintf(stderr, "Error: Compressed input on stdin is not supported. Decompress first:\n");
            fprintf(stderr, "  zcat file.fa.gz | teloscope -o results/\n");
            exit(EXIT_FAILURE);
        }
    }

    // no input
    if (userInput.inSequence.empty() && userInput.pipeType == 'n') {
        fprintf(stderr, "Error: No input file provided. Use -f or pass as positional argument.\n");
        exit(EXIT_FAILURE);
    }

    // step must be <= window
    if (userInput.step > userInput.windowSize) {
        fprintf(stderr, "Error: Step size (%d) cannot be larger than window size (%d).\n",
                userInput.step, userInput.windowSize);
        exit(EXIT_FAILURE);
    } else if (userInput.step < userInput.windowSize) {
        fprintf(stderr, "Sliding with window size (%d) and step size (%d).\n",
                userInput.windowSize, userInput.step);
    }

    // writable check
    if (!userInput.outRoute.empty()) {
        std::string testPath = userInput.outRoute + "/.teloscope_write_test";
        std::ofstream test(testPath);
        if (!test.is_open()) {
            fprintf(stderr, "Error: Output directory '%s' is not writable.\n",
                    userInput.outRoute.c_str());
            exit(EXIT_FAILURE);
        }
        test.close();
        std::filesystem::remove(testPath);
    }

    // default to canonical if -p not provided
    if (!hasInputPatterns) {
        userInput.rawPatterns = {userInput.canonicalFwd, userInput.canonicalRev};
    }
    if (userInput.rawPatterns.empty()) {
        fprintf(stderr, "Warning: No valid patterns supplied via -p. Using canonical: %s, %s\n",
                userInput.canonicalFwd.c_str(), userInput.canonicalRev.c_str());
        userInput.rawPatterns = {userInput.canonicalFwd, userInput.canonicalRev};
    }

    lg.verbose("Input variables assigned");
    userInput.patternInfo = expandPatternsWithOrientation(
        userInput.rawPatterns, userInput.editDistance, userInput.canonicalFwd);
    
    userInput.patterns.clear();
    userInput.patterns.reserve(userInput.patternInfo.size());
    for (const auto& [pattern, isForward] : userInput.patternInfo) {
        userInput.patterns.push_back(pattern);
    }

    fprintf(stderr, "Scanning %zu telomeric variants (includes reverse complements).\n",
            userInput.patterns.size());
    if (userInput.editDistance > 0) {
        fprintf(stderr, "Edit distance enabled: up to %u substitution%s per seed.\n",
                userInput.editDistance,
                userInput.editDistance > 1 ? "s" : "");
    }
    if (userInput.patterns.size() > 500) {
        fprintf(stderr, "Warning: %zu patterns is unusually high and may be slow on large genomes.\n",
                userInput.patterns.size());
        fprintf(stderr, "  Consider fewer IUPAC wildcards or a lower -x value.\n");
    }

    // output summary
    std::string outputSummary;
    auto appendOutput = [&](const char *label) {
        if (!outputSummary.empty()) {
            outputSummary += ", ";
        }
        outputSummary += label;
    };
    if (userInput.outGC) appendOutput("GC windows");
    if (userInput.outWinRepeats) appendOutput("repeat density");
    if (userInput.outEntropy) appendOutput("Shannon entropy");
    if (userInput.outMatches) appendOutput("genome-wide matches");
    if (userInput.outITS) appendOutput("ITS blocks");
    if (userInput.outPlotReport) appendOutput("plot report");
    if (!outputSummary.empty()) {
        fprintf(stderr, "Outputs: %s.\n", outputSummary.c_str());
    }

    // command echo
    if (cmd_flag) {
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }

    // Start processing threads and load inputs
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

    if (userInput.outPlotReport) {
        std::string scriptPath = findReportScript(argv[0]);
        if (scriptPath.empty()) {
            fprintf(stderr, "Warning: Could not locate teloscope_report.py.\n");
            fprintf(stderr, "  Ensure the scripts/ directory is present alongside the teloscope binary.\n");
        } else {
            std::string pdfPath = userInput.outRoute + "/" + userInput.inSequenceName + "_plot_report.pdf";
            fprintf(stderr, "Generating plot report: %s\n", pdfPath.c_str());

            std::string cmd = "python3 \"" + scriptPath + "\" \""
                            + userInput.outRoute + "\" -o \"" + pdfPath + "\"";
            int ret = system(cmd.c_str());
            if (ret != 0) {
                fprintf(stderr, "Warning: Report generation failed.\n");
                fprintf(stderr, "  Ensure Python 3 is installed with: matplotlib, numpy, pandas\n");
                fprintf(stderr, "  Or generate manually: python3 scripts/teloscope_report.py %s\n",
                        userInput.outRoute.c_str());
            }
        }
    }

    exit(EXIT_SUCCESS);

}
