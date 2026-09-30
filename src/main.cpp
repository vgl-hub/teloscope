#include "main.h"
#include "bam.h"
#include "read-filter.h"
#include "teloscope.h"
#include <input.h>
#include <iostream>
#include <zlib.h>

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

#ifdef _WIN32
#include <fcntl.h>
#endif

std::string version = "0.1.6";

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

// backs std::cin with fread(stdin) so detection and every reader avoid the default per-character stdio sync cost
class StdinBuffer : public std::streambuf {
public:
    StdinBuffer() { setg(buffer, buffer, buffer); }
protected:
    int_type underflow() override {
        if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
        const size_t got = fread(buffer, 1, sizeof(buffer), stdin);
        if (got == 0) return traits_type::eof();
        setg(buffer, buffer, buffer + got);
        return traits_type::to_int_type(*gptr());
    }
private:
    char buffer[1 << 16]; // 64 KiB
};

// an assembly must start like FASTA ('>') or GFA ('#', or a record letter and a tab)
static bool looksLikeAssembly(unsigned char first, unsigned char second) {
    return first == '>' || first == '#' ||
        (std::string("HSLCPWJEFGOU").find(static_cast<char>(first)) != std::string::npos && second == '\t');
}

// FASTQ or BAM input: the telomeric reads, the read telomere BED and the report, all under outRoute
static void runReads(Input &in) {
    const bool bam = userInput.readInput == ReadInput::bam;
    const std::string base = userInput.outRoute + "/" + userInput.inSequenceName;
    const std::string bedPath = base + "_terminal_telomeres.bed", reportPath = base + "_report.tsv";
    const std::string subsetPath = bam
        ? userInput.outRoute + "/" + std::filesystem::path(userInput.inSequenceName).stem().string() + "_telomeric.bam"
        : base + "_telomeric.fastq";

    std::ofstream subset(subsetPath, std::ios::binary), bed(bedPath), report(reportPath);
    auto removeOutputs = [&]() {
        subset.close(); bed.close(); report.close(); // Windows cannot remove an open file
        for (const std::string &path : {subsetPath, bedPath, reportPath}) {
            std::error_code removeError;
            std::filesystem::remove(path, removeError);
        }
    };
    if (!subset.is_open() || !bed.is_open() || !report.is_open()) {
        fprintf(stderr, "Error: cannot write telomeric records to '%s'.\n", userInput.outRoute.c_str());
        removeOutputs();
        threadPool.join();
        exit(EXIT_FAILURE);
    }

    ReadTlStats stats;
    try {
        if (bam) readBamReads(userInput, subset, bed, stats);
        else in.readFastqReads(subset, bed, stats);
        if (!bed.flush().good()) throw std::runtime_error("failed while writing the read telomere BED");
        writeReadTlReport(report, makeReadTlInput(userInput), stats);
        if (!report.flush().good()) throw std::runtime_error("failed while writing the report");
    } catch (const std::exception &error) {
        fprintf(stderr, "Error: %s.\n", error.what());
        if (bam && userInput.inSequence.empty())
            fprintf(stderr, "  Compressed FASTQ or FASTA on stdin or a pipe is not supported: pass the file, or decompress it first (zcat reads.fq.gz | teloscope).\n");
        threadPool.join();
        removeOutputs();
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "Wrote %s, %s and %s.\n", subsetPath.c_str(), bedPath.c_str(), reportPath.c_str());
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
    bool fullScanRequested = false;
    
    if (argc == 1 && !isPipe) { // case: with no arguments and no pipe

        printf("teloscope input.[fa|fa.gz|gfa|fq|fq.gz|bam] [options]\nUse -h for additional help.\n");
        exit(0);

    }

    auto setInputFile = [&](const char* path) {
        ifFileExists(path);
        userInput.inSequence = path;
        std::error_code pathError; // pipes and fd paths have no canonical form
        const std::filesystem::path real = std::filesystem::canonical(userInput.inSequence, pathError);
        const std::filesystem::path resolved = pathError ? std::filesystem::path(path) : real;
        userInput.inSequence       = resolved.string();
        userInput.inSequencePrefix = resolved.parent_path().string();
        userInput.inSequenceName   = resolved.filename().string();
        if (!userInput.outRouteSet) // an explicit -o is not overridden by the input directory
            userInput.outRoute = userInput.inSequencePrefix;
    };

    auto addBedFilterFile = [&](const char* path, std::vector<std::string>& files,
                                const char* optionName) {
        const std::filesystem::path selectorPath(path ? path : "");
        std::error_code error;
        if (selectorPath.empty() || !std::filesystem::exists(selectorPath, error) || error) {
            fprintf(stderr, "Error: %s file does not exist: '%s'.\n", optionName,
                    path ? path : "");
            exit(EXIT_FAILURE);
        }
        if (!std::filesystem::is_regular_file(selectorPath, error) || error) {
            fprintf(stderr, "Error: %s file '%s' is not a regular file.\n", optionName, path);
            exit(EXIT_FAILURE);
        }
        const std::filesystem::path resolved = std::filesystem::canonical(selectorPath, error);
        if (error) {
            fprintf(stderr, "Error: Could not resolve %s file '%s': %s.\n", optionName,
                    path, error.message().c_str());
            exit(EXIT_FAILURE);
        }
        files.push_back(resolved.string());
        userInput.sequenceFilterActive = true;
    };

    auto addPrefixFilters = [&](const char* value, std::vector<std::string>& prefixes,
                                const char* optionName) {
        std::istringstream stream(value ? value : "");
        std::string prefix;
        bool found = false;
        while (std::getline(stream, prefix, ',')) {
            const size_t first = prefix.find_first_not_of(" \t\r\n");
            const size_t last = prefix.find_last_not_of(" \t\r\n");
            if (first == std::string::npos) {
                fprintf(stderr, "Error: %s contains an empty prefix.\n", optionName);
                exit(EXIT_FAILURE);
            }
            prefix = prefix.substr(first, last - first + 1);
            prefixes.push_back(prefix);
            found = true;
        }
        if (!found || (value && value[0] != '\0' && value[strlen(value) - 1] == ',')) {
            fprintf(stderr, "Error: %s contains an empty prefix.\n", optionName);
            exit(EXIT_FAILURE);
        }
        userInput.sequenceFilterActive = true;
    };

    // true when the whole value is one number: 5e4, 10kb or 1.5 for an integer stop at the letter
    auto wholeLong = [](const char* value, long& v) -> bool {
        size_t used = 0;
        try { v = std::stol(value, &used); } catch (const std::exception&) { return false; }
        return value[used] == '\0';
    };
    auto wholeFloat = [](const char* value, float& v) -> bool {
        size_t used = 0;
        try { v = std::stof(value, &used); } catch (const std::exception&) { return false; }
        return value[used] == '\0';
    };

    auto parsePositive = [&](const char* value, const char* optionName) -> uint32_t {
        long v = 0;
        if (!wholeLong(value, v) || v > INT32_MAX) {
            fprintf(stderr, "Error: Invalid value '%s' for %s. Must be a whole number up to 2147483647.\n", value, optionName);
            exit(EXIT_FAILURE);
        }
        if (v <= 0) {
            fprintf(stderr, "Error: %s must be > 0.\n", optionName);
            exit(EXIT_FAILURE);
        }
        return static_cast<uint32_t>(v);
    };

    auto parseLabelThreshold = [&](const char* value) -> float {
        float v = 0;
        if (!wholeFloat(value, v)) {
            fprintf(stderr, "Error: Invalid value '%s' for --label-threshold. Must be a number.\n", value);
            exit(EXIT_FAILURE);
        }
        if (v <= 0.5f || v > 1.0f) {
            fprintf(stderr, "Error: --label-threshold must be in the range (0.5,1].\n");
            exit(EXIT_FAILURE);
        }
        return v;
    };

    static struct option long_options[] = { // struct mapping long options
        {"input-sequence", required_argument, 0, 'f'},
        {"output", required_argument, 0, 'o'},
        {"include-bed", required_argument, 0, 0},
        {"exclude-bed", required_argument, 0, 0},
        {"include-prefix", required_argument, 0, 0},
        {"exclude-prefix", required_argument, 0, 0},
        {"chr-only", no_argument, 0, 0},
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
        {"terminal-tolerance", required_argument, 0, 0},
        {"label-threshold", required_argument, 0, 0},
        {"min-block-counts", required_argument, 0, 0},

        {"out-fasta", no_argument, 0, 'a'},
        {"out-win-repeats", no_argument, 0, 'r'},
        {"out-gc", no_argument, 0, 'g'},
        {"out-entropy", no_argument, 0, 'e'},
        {"out-matches", no_argument, 0, 'm'},
        {"out-its", no_argument, 0, 'i'},
        {"ultra-fast", no_argument, 0, 'u'},
        {"manual-curation", no_argument, 0, 'n'},
        {"plot-report", no_argument, 0, 0},
        {"fastq-subset", no_argument, 0, 0},
        {"bam-subset", no_argument, 0, 0},
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
                if (optopt == 0 && optind > 0 && optind <= argc &&
                    strncmp(argv[optind - 1], "--", 2) == 0) {
                    fprintf(stderr, "Error: Option %s is missing a required argument\n",
                            argv[optind - 1]);
                } else {
                    fprintf(stderr, "Error: Option -%c is missing a required argument\n", optopt);
                }
                return EXIT_FAILURE;
            default: // unknown or ambiguous option, before any input is read
                if (optopt == 0 && optind > 0 && optind <= argc &&
                    strncmp(argv[optind - 1], "--", 2) == 0) {
                    fprintf(stderr, "Error: Unknown or ambiguous option %s\n", argv[optind - 1]);
                } else {
                    fprintf(stderr, "Error: Unknown option -%c\n", optopt);
                }
                return EXIT_FAILURE;


            case 0: // long options without short options
                if (strcmp(long_options[option_index].name, "plot-report") == 0)
                    userInput.outPlotReport = true;
                else if (strcmp(long_options[option_index].name, "fastq-subset") == 0 ||
                         strcmp(long_options[option_index].name, "bam-subset") == 0) {
                    fprintf(stderr, "Error: --fastq-subset and --bam-subset were removed: FASTQ or BAM input is detected and always writes the telomeric reads, the read telomere BED and the report.\n");
                    exit(EXIT_FAILURE);
                }
                else if (strcmp(long_options[option_index].name, "include-bed") == 0)
                    addBedFilterFile(optarg, userInput.includeBedFiles, "--include-bed");
                else if (strcmp(long_options[option_index].name, "exclude-bed") == 0)
                    addBedFilterFile(optarg, userInput.excludeBedFiles, "--exclude-bed");
                else if (strcmp(long_options[option_index].name, "include-prefix") == 0)
                    addPrefixFilters(optarg, userInput.includePrefixes, "--include-prefix");
                else if (strcmp(long_options[option_index].name, "exclude-prefix") == 0)
                    addPrefixFilters(optarg, userInput.excludePrefixes, "--exclude-prefix");
                else if (strcmp(long_options[option_index].name, "chr-only") == 0) {
                    userInput.chrOnly = true;
                    userInput.sequenceFilterActive = true;
                }
                else if (strcmp(long_options[option_index].name, "terminal-tolerance") == 0) {
                    userInput.terminalTolerance = parsePositive(optarg, "--terminal-tolerance");
                    userInput.terminalToleranceSet = true;
                }
                else if (strcmp(long_options[option_index].name, "label-threshold") == 0)
                    userInput.labelThreshold = parseLabelThreshold(optarg);
                else if (strcmp(long_options[option_index].name, "min-block-counts") == 0)
                    userInput.minBlockCounts = parsePositive(optarg, "--min-block-counts");
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
                    userInput.outRouteSet = true;

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


            case 'w':
                userInput.windowSize = parsePositive(optarg, "-w/--window");
                break;


            case 's':
                userInput.step = parsePositive(optarg, "-s/--step");
                break;


            case 't':
                userInput.terminalLimit = parsePositive(optarg, "-t/--terminal-limit");
                userInput.terminalLimitSet = true;
                break;


            case 'k': // max match distance
                userInput.maxMatchDist = parsePositive(optarg, "-k/--max-match-distance");
                break;


            case 'l':
                userInput.minBlockLen = parsePositive(optarg, "-l/--min-block-length");
                break;


            case 'd':
                userInput.maxBlockDist = parsePositive(optarg, "-d/--max-block-distance");
                break;


            case 'y': { // min block density
                float v = 0;
                if (!wholeFloat(optarg, v)) {
                    fprintf(stderr, "Error: Invalid min block density '%s'. Must be a number in (0,1].\n", optarg);
                    exit(EXIT_FAILURE);
                }
                if (v <= 0.0f || v > 1.0f) {
                    fprintf(stderr, "Error: Min block density (-y/--min-block-density) must be in the range (0,1].\n");
                    exit(EXIT_FAILURE);
                }
                userInput.minBlockDensity = v;
                break;
            }


            case 'x': { // edit distance
                long v = 0;
                if (!wholeLong(optarg, v)) {
                    fprintf(stderr, "Error: Invalid edit distance '%s'. Must be a number [0,2].\n", optarg);
                    exit(EXIT_FAILURE);
                }
                if (v < 0 || v > 2) {
                    fprintf(stderr, "Error: Edit distance (-x/--edit-distance) must be in the range [0,2].\n");
                    exit(EXIT_FAILURE);
                }
                userInput.editDistance = static_cast<uint8_t>(v);
                break;
            }


            case 'a':
                userInput.outFasta = true;
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
                fullScanRequested = true;
                userInput.ultraFastMode = false;
                break;


            case 'u': {
                if (userInput.outWinRepeats || userInput.outGC ||
                    userInput.outEntropy   || fullScanRequested ||
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


            case 'n': // manual curation mode: also report contig-internal telomeres
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
                printf("teloscope input.[fa|fa.gz|gfa|fq|fq.gz|bam] [options]\n");
                printf("teloscope -f input.[fa|fa.gz|gfa|fq|fq.gz|bam] [options]\n");
                printf("FASTQ or BAM input is detected and writes the telomeric reads, the read telomere BED and the report (estimates; see docs).\n");
                printf("\nRequired Parameters:\n");
                printf("\t'-f'\t--input-sequence\tInput FASTA, GFA, FASTQ, or BAM file (or pass as first positional argument).\n");
                printf("\t'-o'\t--output\tSet output route. [Default: Input path]\n");
                printf("\t'-c'\t--canonical\tSet canonical pattern. [Default: TTAGGG]\n");
                printf("\t'-p'\t--patterns\tSet patterns to explore, separate them by commas [Default: TTAGGG]\n");
                printf("\t'-j'\t--threads\tSet maximum number of threads. [Default: max. available]\n");
                printf("\t'-t'\t--terminal-limit\tSet how far in from each end to look; also the read tile size. [Default: 50000 assembly, 2000 reads]\n");
                printf("\t'-k'\t--max-match-distance\tSet maximum distance for merging matches. [Default: 50]\n");
                printf("\t'-d'\t--max-block-distance\tSet maximum non-telomeric stretch inside a telomere. [Default: 1000]\n");
                printf("\t'-l'\t--min-block-length\tSet minimum block length. [Default: 300]\n");
                printf("\t'-y'\t--min-block-density\tSet minimum block density. [Default: 0.5]\n");
                printf("\t'-x'\t--edit-distance\tSet edit distance for pattern matching (0-2). [Default: 1]\n");
                printf("\t\t--terminal-tolerance\tSet how far from an end a telomere may start. [Default: 3000 assembly, 300 reads]\n");
                printf("\t\t--label-threshold\tSet forward-strand fraction for the p/q label. [Default: 0.667]\n");
                printf("\t\t--min-block-counts\tSet minimum canonical matches per block. [Default: 2]\n");

                printf("\nOptional Parameters:\n");
                printf("\t'-w'\t--window\tSet sliding window size. [Default: 1000]\n");
                printf("\t'-s'\t--step\tSet sliding window step. [Default: 1000 (non-overlapping)]\n");
                printf("\t\t--include-bed FILE\tAnalyze whole records whose IDs occur in BED/list column 1. Repeatable. [Default: unset]\n");
                printf("\t\t--exclude-bed FILE\tExclude whole records whose IDs occur in BED/list column 1. Repeatable. [Default: unset]\n");
                printf("\t\t--include-prefix LIST\tAnalyze IDs with a literal, comma-separated prefix. Repeatable. [Default: unset]\n");
                printf("\t\t--exclude-prefix LIST\tExclude IDs with a literal, comma-separated prefix. Repeatable. [Default: unset]\n");
                printf("\t\t--chr-only\tKeep only records named like the longest one (chromosome convention). [Default: false]\n");
                printf("\t\tRecord filters are off by default. Matching is case-sensitive; includes form a union and exclusions apply last. FASTA uses the first ID token; BED never crops records.\n");
                printf("\t'-r'\t--out-win-repeats\tOutput per-window repeat density, canonical ratio, and strand ratio. [Default: false]\n");
                printf("\t'-g'\t--out-gc\tOutput GC content for each window. [Default: false]\n");
                printf("\t'-e'\t--out-entropy\tOutput Shannon entropy for each window. [Default: false]\n");
                printf("\t'-m'\t--out-matches\tOutput all canonical and terminal non-canonical matches. [Default: false]\n");
                printf("\t'-a'\t--out-fasta\tOutput terminal telomere sequences as FASTA. [Default: false]\n");
                printf("\t'-i'\t--out-its\tScan whole sequences for interstitial telomeres. [Default: false]\n");
                printf("\t'-u'\t--ultra-fast\tUltra-fast mode. Only scans terminal telomeres at scaffold ends. [Default: true]\n");
                printf("\t'-n'\t--manual-curation\tAlso report telomeres at contig ends. [Default: scaffold only]\n");
                printf("\t\t--plot-report\tGenerate terminal and ITS PDF reports after analysis (requires Python 3 + matplotlib). [Default: false]\n");

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
        if (!userInput.outRouteSet)
            userInput.outRoute = ".";
    }

    // no input
    if (userInput.inSequence.empty() && userInput.pipeType == 'n') {
        fprintf(stderr, "Error: No input file provided. Use -f or pass as positional argument.\n");
        exit(EXIT_FAILURE);
    }

    // reads are detected by content: the BAM magic, a FASTQ header, or a FASTA/GFA start; anything else is refused
    std::error_code fileError;
    if (std::filesystem::is_directory(userInput.inSequence, fileError)) {
        fprintf(stderr, "Error: input '%s' is a directory.\n", userInput.inSequence.c_str());
        exit(EXIT_FAILURE);
    }

    // a named path that is not a regular file (a FIFO, /dev/fd/N, /dev/stdin, a character device) reads exactly like stdin
    std::string namedPipeInput;
    if (!userInput.inSequence.empty() && !std::filesystem::is_regular_file(userInput.inSequence, fileError)) {
        namedPipeInput = userInput.inSequence;
        if (!freopen(namedPipeInput.c_str(), "rb", stdin)) {
            fprintf(stderr, "Error: cannot open input '%s'.\n", namedPipeInput.c_str());
            exit(EXIT_FAILURE);
        }
        userInput.inSequence.clear();
        userInput.pipeType = 'f';
        if (!userInput.outRouteSet)
            userInput.outRoute = ".";
    }

    if (userInput.inSequence.empty()) {
        static StdinBuffer stdinBuffer;
        std::cin.rdbuf(&stdinBuffer); // gfalibs, bam.cpp and the filtered FASTA loader all read std::cin
#ifdef _WIN32
        _setmode(_fileno(stdin), _O_BINARY); // detection buffers bytes before bam.cpp would otherwise set this
#endif
        const int first = std::cin.peek(); // gzip on stdin/a pipe can only be a BAM; BgzfReader checks the magic
        if (first == EOF) {
            if (namedPipeInput.empty())
                fprintf(stderr, "Error: input on stdin is empty.\n");
            else
                fprintf(stderr, "Error: input '%s' is empty.\n", namedPipeInput.c_str());
            exit(EXIT_FAILURE);
        }
        if (first == '@') userInput.readInput = ReadInput::fastq;
        else if (first == 0x1f) userInput.readInput = ReadInput::bam;
        else {
            std::cin.get(); // both bytes sit in stdinBuffer, so unget restores the first
            const int second = std::cin.peek();
            std::cin.unget();
            if (!looksLikeAssembly(first, second)) {
                if (namedPipeInput.empty())
                    fprintf(stderr, "Error: input on stdin is not FASTA, GFA, FASTQ or BAM.\n");
                else
                    fprintf(stderr, "Error: input '%s' is not FASTA, GFA, FASTQ or BAM.\n", namedPipeInput.c_str());
                exit(EXIT_FAILURE);
            }
        }
    } else {
        char magic[4] = {0, 0, 0, 0}; // gzread inflates gzip/BGZF and passes plain files through
        if (gzFile file = gzopen(userInput.inSequence.c_str(), "rb")) {
            const int got = gzread(file, magic, sizeof(magic));
            int status = Z_OK;
            gzerror(file, &status);
            gzclose(file);
            if (status != Z_OK && status != Z_STREAM_END) {
                fprintf(stderr, "Error: cannot inflate compressed input '%s'.\n", userInput.inSequence.c_str());
                exit(EXIT_FAILURE);
            }
            if (got == 0) {
                fprintf(stderr, "Error: input '%s' is empty.\n", userInput.inSequence.c_str());
                exit(EXIT_FAILURE);
            }
        } else {
            fprintf(stderr, "Error: cannot open input '%s'.\n", userInput.inSequence.c_str());
            exit(EXIT_FAILURE);
        }
        if (memcmp(magic, "BAM\1", 4) == 0) userInput.readInput = ReadInput::bam;
        else if (magic[0] == '@') userInput.readInput = ReadInput::fastq;
        else if (!looksLikeAssembly(magic[0], magic[1])) {
            fprintf(stderr, "Error: input '%s' is not FASTA, GFA, FASTQ or BAM.\n", userInput.inSequence.c_str());
            exit(EXIT_FAILURE);
        }
    }

    if (userInput.sequenceFilterActive && userInput.readInput != ReadInput::none) {
        fprintf(stderr, "Error: --include-bed/--exclude-bed/--include-prefix/--exclude-prefix/--chr-only "
                        "filter assembly records and cannot be used with FASTQ or BAM input.\n");
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
    // only the assembly FASTA path scans in windows, so only it can outgrow one
    const bool usesWindows = !userInput.ultraFastMode && userInput.readInput == ReadInput::none &&
                             !isGfaAssemblyPath(userInput.inSequence);
    for (const auto& [pattern, isForward] : userInput.patternInfo) {
        userInput.patterns.push_back(pattern);
        if (pattern.size() > 255) { // a match records its size in eight bits
            fprintf(stderr, "Error: Pattern '%s' is longer than 255 bases.\n", pattern.c_str());
            exit(EXIT_FAILURE);
        }
        if (usesWindows && pattern.size() > userInput.windowSize) {
            fprintf(stderr, "Error: Window size (%u) is smaller than pattern '%s'.\n",
                    userInput.windowSize, pattern.c_str());
            exit(EXIT_FAILURE);
        }
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
    if (userInput.outFasta) appendOutput("telomere FASTA");
    if (fullScanRequested) appendOutput("full scan");
    if (userInput.outPlotReport) appendOutput("plot report");
    if (!outputSummary.empty()) {
        fprintf(stderr, "Outputs: %s.\n", outputSummary.c_str());
    }
    if (userInput.readInput != ReadInput::none &&
        (userInput.outFasta || userInput.outWinRepeats || userInput.outGC ||
         userInput.outEntropy || userInput.outMatches || fullScanRequested ||
         userInput.outPlotReport || userInput.manualCuration)) {
        fprintf(stderr, "Warning: assembly output flags are ignored for FASTQ or BAM input.\n");
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

    if (userInput.readInput != ReadInput::none) {
        runReads(in);
        threadPool.join();
        exit(EXIT_SUCCESS);
    }

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
            std::cout.flush();

            // the run's file stem names this input's files, so other runs in the directory are left alone
            std::string cmd = "python3 \"" + scriptPath + "\" \"" + userInput.outRoute + "/" + userInput.inSequenceName + "\"";
            int ret = system(cmd.c_str());
            if (ret != 0) {
                fprintf(stderr, "Warning: Report generation failed.\n");
                fprintf(stderr, "  Ensure Python 3 is installed with: matplotlib, numpy, pandas\n");
                fprintf(stderr, "  Or generate manually: python3 scripts/teloscope_report.py %s/%s\n",
                        userInput.outRoute.c_str(), userInput.inSequenceName.c_str());
            }
        }
    }

    exit(EXIT_SUCCESS);

}
