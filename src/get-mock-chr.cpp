#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>

// ── structs ──────────────────────────────────────────────────────────────────

struct SimParams {
    uint32_t n        = 1000;   // sequences per run
    uint32_t R        = 2000;   // canonical repeats per tip
    uint32_t T        = 100;    // TVR repeats per tip
    uint32_t L        = 25000;  // central segment length
    double   rate     = 0.0;    // per-base mutation rate
    uint64_t seed     = 42;
    std::string outDir = "testFiles/simulate";
};

struct GroundTruth {
    std::string name;
    uint64_t len;
    uint64_t pTeloStart, pTeloEnd;
    uint64_t pTVRStart,  pTVREnd;
    uint64_t qTVRStart,  qTVREnd;
    uint64_t qTeloStart, qTeloEnd;
};

struct BEDEntry {
    std::string chrom;
    uint64_t start, end;
    char label; // p, q, b
};

// ── helpers ──────────────────────────────────────────────────────────────────

static const char BASES[] = "ACGT";

static char randomBase(std::mt19937_64 &rng) {
    return BASES[std::uniform_int_distribution<int>(0, 3)(rng)];
}

static char mutateBase(char orig, std::mt19937_64 &rng) {
    char b;
    do { b = randomBase(rng); } while (b == orig);
    return b;
}

static void mkdirP(const std::string &path) {
    std::filesystem::create_directories(path);
}

static void writeFastaEntry(std::ofstream &out, const std::string &name, const std::string &seq) {
    out << '>' << name << '\n';
    for (size_t i = 0; i < seq.size(); i += 80)
        out << seq.substr(i, 80) << '\n';
}

// ── generate ─────────────────────────────────────────────────────────────────

static void generate(const SimParams &p) {
    // Seed: combine user seed with hash of rate for per-rate determinism
    std::hash<double> dh;
    uint64_t combinedSeed = p.seed + dh(p.rate);
    std::mt19937_64 rng(combinedSeed);

    mkdirP(p.outDir);

    std::string faPath  = p.outDir + "/sequences.fa";
    std::string tsvPath = p.outDir + "/ground_truth.tsv";

    std::ofstream fa(faPath);
    std::ofstream tsv(tsvPath);
    if (!fa || !tsv) {
        std::cerr << "Error: cannot open output files in " << p.outDir << '\n';
        exit(1);
    }

    // TSV header
    tsv << "# seed=" << p.seed
        << " n=" << p.n
        << " R=" << p.R
        << " T=" << p.T
        << " L=" << p.L
        << " rate=" << p.rate << '\n';
    tsv << "name\tlen\tpTeloStart\tpTeloEnd\tpTVRStart\tpTVREnd"
        << "\tqTVRStart\tqTVREnd\tqTeloStart\tqTeloEnd\n";

    const std::string pCanonUnit = "CCCTAA"; // forward pattern → p-arm
    const std::string qCanonUnit = "TTAGGG"; // reverse pattern → q-arm
    const size_t unitLen = 6;

    std::uniform_real_distribution<double> mutDist(0.0, 1.0);

    for (uint32_t i = 0; i < p.n; ++i) {
        std::string seq;
        // Pre-calculate sizes
        uint64_t pCanonLen = p.R * unitLen;
        uint64_t pTVRLen   = p.T * unitLen;
        uint64_t qTVRLen   = p.T * unitLen;
        uint64_t qCanonLen = p.R * unitLen;
        uint64_t totalLen  = pCanonLen + pTVRLen + p.L + qTVRLen + qCanonLen;
        seq.reserve(totalLen);

        // p-arm canonical: R × CCCTAA
        for (uint32_t r = 0; r < p.R; ++r)
            seq += pCanonUnit;

        // p-arm TVR: T × HD1(CCCTAA) — one random substitution per repeat
        for (uint32_t t = 0; t < p.T; ++t) {
            std::string unit = pCanonUnit;
            int pos = std::uniform_int_distribution<int>(0, unitLen - 1)(rng);
            unit[pos] = mutateBase(unit[pos], rng);
            seq += unit;
        }

        // Central: random ACGT
        for (uint32_t c = 0; c < p.L; ++c)
            seq += randomBase(rng);

        // q-arm TVR: T × HD1(TTAGGG)
        for (uint32_t t = 0; t < p.T; ++t) {
            std::string unit = qCanonUnit;
            int pos = std::uniform_int_distribution<int>(0, unitLen - 1)(rng);
            unit[pos] = mutateBase(unit[pos], rng);
            seq += unit;
        }

        // q-arm canonical: R × TTAGGG
        for (uint32_t r = 0; r < p.R; ++r)
            seq += qCanonUnit;

        // Apply global mutations at rate r
        if (p.rate > 0.0) {
            for (size_t j = 0; j < seq.size(); ++j) {
                if (mutDist(rng) < p.rate)
                    seq[j] = mutateBase(seq[j], rng);
            }
        }

        // Ground truth coordinates
        uint64_t pTeloStart = 0;
        uint64_t pTeloEnd   = pCanonLen;
        uint64_t pTVRStart  = pCanonLen;
        uint64_t pTVREnd    = pCanonLen + pTVRLen;
        uint64_t qTVRStart  = pCanonLen + pTVRLen + p.L;
        uint64_t qTVREnd    = qTVRStart + qTVRLen;
        uint64_t qTeloStart = qTVREnd;
        uint64_t qTeloEnd   = totalLen;

        std::string name = "sim_" + std::to_string(i + 1);

        writeFastaEntry(fa, name, seq);

        tsv << name << '\t' << totalLen
            << '\t' << pTeloStart << '\t' << pTeloEnd
            << '\t' << pTVRStart  << '\t' << pTVREnd
            << '\t' << qTVRStart  << '\t' << qTVREnd
            << '\t' << qTeloStart << '\t' << qTeloEnd << '\n';
    }

    std::cerr << "Generated " << p.n << " sequences in " << p.outDir << '\n';
}

// ── evaluate ─────────────────────────────────────────────────────────────────

static std::vector<GroundTruth> readGroundTruth(const std::string &path) {
    std::vector<GroundTruth> gts;
    std::ifstream in(path);
    if (!in) { std::cerr << "Error: cannot open " << path << '\n'; exit(1); }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.substr(0, 4) == "name") continue;
        std::istringstream ss(line);
        GroundTruth g;
        ss >> g.name >> g.len
           >> g.pTeloStart >> g.pTeloEnd
           >> g.pTVRStart  >> g.pTVREnd
           >> g.qTVRStart  >> g.qTVREnd
           >> g.qTeloStart >> g.qTeloEnd;
        gts.push_back(g);
    }
    return gts;
}

static std::vector<BEDEntry> readBED(const std::string &path) {
    std::vector<BEDEntry> entries;
    std::ifstream in(path);
    if (!in) { std::cerr << "Error: cannot open " << path << '\n'; exit(1); }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        BEDEntry e;
        std::string labelStr;
        uint32_t blockLen;
        // BED format: chrom start end blockLen label ...
        ss >> e.chrom >> e.start >> e.end >> blockLen >> labelStr;
        e.label = labelStr.empty() ? '?' : labelStr[0];
        entries.push_back(e);
    }
    return entries;
}

static bool overlaps(uint64_t s1, uint64_t e1, uint64_t s2, uint64_t e2) {
    return s1 < e2 && s2 < e1;
}

static void evaluate(const std::string &gtPath, const std::string &bedPath) {
    auto gts = readGroundTruth(gtPath);
    auto beds = readBED(bedPath);

    uint64_t totalTips = gts.size() * 2; // p + q for each contig
    uint64_t detected = 0;
    double sumBias = 0.0;
    double sumAbsErr = 0.0;
    uint64_t measuredTips = 0;
    uint64_t tvrInclusions = 0;
    uint64_t tvrMeasured = 0;
    uint64_t fpBlocks = 0;

    // Index BED entries by chrom
    std::vector<bool> bedUsed(beds.size(), false);

    for (const auto &g : gts) {
        bool pDetected = false, qDetected = false;

        for (size_t bi = 0; bi < beds.size(); ++bi) {
            const auto &b = beds[bi];
            if (b.chrom != g.name) continue;

            // Check p-arm: label p or b, overlapping canonical region
            if ((b.label == 'p' || b.label == 'b') &&
                overlaps(b.start, b.end, g.pTeloStart, g.pTeloEnd)) {
                pDetected = true;
                bedUsed[bi] = true;

                // Tract-length error: estimated vs true canonical length
                uint64_t estLen = b.end - b.start;
                uint64_t trueLen = g.pTeloEnd - g.pTeloStart;
                sumBias += (double)estLen - (double)trueLen;
                sumAbsErr += std::abs((double)estLen - (double)trueLen);
                measuredTips++;

                // TVR inclusion: does BED block extend into TVR?
                if (overlaps(b.start, b.end, g.pTVRStart, g.pTVREnd))
                    tvrInclusions++;
                tvrMeasured++;
            }

            // Check q-arm: label q or b, overlapping canonical region
            if ((b.label == 'q' || b.label == 'b') &&
                overlaps(b.start, b.end, g.qTeloStart, g.qTeloEnd)) {
                qDetected = true;
                bedUsed[bi] = true;

                uint64_t estLen = b.end - b.start;
                uint64_t trueLen = g.qTeloEnd - g.qTeloStart;
                sumBias += (double)estLen - (double)trueLen;
                sumAbsErr += std::abs((double)estLen - (double)trueLen);
                measuredTips++;

                if (overlaps(b.start, b.end, g.qTVRStart, g.qTVREnd))
                    tvrInclusions++;
                tvrMeasured++;
            }
        }

        if (pDetected) detected++;
        if (qDetected) detected++;
    }

    // Count false positives: BED entries not matched to any truth region
    for (size_t bi = 0; bi < beds.size(); ++bi) {
        if (bedUsed[bi]) continue;
        // Check if this entry overlaps ANY true telomere region
        bool matchesAny = false;
        for (const auto &g : gts) {
            if (beds[bi].chrom != g.name) continue;
            if (overlaps(beds[bi].start, beds[bi].end, g.pTeloStart, g.pTVREnd) ||
                overlaps(beds[bi].start, beds[bi].end, g.qTVRStart, g.qTeloEnd)) {
                matchesAny = true;
                break;
            }
        }
        if (!matchesAny) fpBlocks++;
    }

    double sensitivity = totalTips > 0 ? (double)detected / totalTips : 0.0;
    double meanBias    = measuredTips > 0 ? sumBias / measuredTips : 0.0;
    double meanAbsErr  = measuredTips > 0 ? sumAbsErr / measuredTips : 0.0;
    double tvrRate     = tvrMeasured > 0 ? (double)tvrInclusions / tvrMeasured : 0.0;

    // Output: total_tips detected sensitivity mean_bias_bp mean_abs_err_bp tvr_rate fp_blocks
    std::cout << totalTips << '\t'
              << detected << '\t'
              << sensitivity << '\t'
              << meanBias << '\t'
              << meanAbsErr << '\t'
              << tvrRate << '\t'
              << fpBlocks;
}

// ── usage ────────────────────────────────────────────────────────────────────

static void printUsage() {
    std::cerr <<
        "Usage:\n"
        "  Generate: teloscope-simulate [options]\n"
        "  Evaluate: teloscope-simulate --evaluate -g <truth.tsv> -b <telomeres.bed>\n"
        "\nGenerate options:\n"
        "  -n INT    Sequences per run       [1000]\n"
        "  -R INT    Canonical repeats/tip   [2000]\n"
        "  -T INT    TVR repeats/tip         [100]\n"
        "  -L INT    Central segment length  [25000]\n"
        "  -r FLOAT  Mutation rate (per bp)  [0.0]\n"
        "  -s INT    Random seed             [42]\n"
        "  -o DIR    Output directory        [testFiles/simulate]\n"
        "\nEvaluate options:\n"
        "  --evaluate   Enable evaluate mode\n"
        "  -g FILE      Ground truth TSV\n"
        "  -b FILE      Terminal telomeres BED\n"
        "  -h           Print this help\n";
}

// ── main ─────────────────────────────────────────────────────────────────────

int main(int argc, char **argv) {
    SimParams params;
    bool evalMode = false;
    std::string gtFile, bedFile;

    static struct option longOpts[] = {
        {"evaluate", no_argument,       nullptr, 'E'},
        {"help",     no_argument,       nullptr, 'h'},
        {nullptr,    0,                 nullptr,  0 }
    };

    int c;
    while ((c = getopt_long(argc, argv, "n:R:T:L:r:s:o:g:b:h", longOpts, nullptr)) != -1) {
        switch (c) {
            case 'n': params.n    = std::stoul(optarg); break;
            case 'R': params.R    = std::stoul(optarg); break;
            case 'T': params.T    = std::stoul(optarg); break;
            case 'L': params.L    = std::stoul(optarg); break;
            case 'r': params.rate = std::stod(optarg);  break;
            case 's': params.seed = std::stoull(optarg); break;
            case 'o': params.outDir = optarg;            break;
            case 'g': gtFile  = optarg;                  break;
            case 'b': bedFile = optarg;                  break;
            case 'E': evalMode = true;                   break;
            case 'h': printUsage(); return 0;
            default:  printUsage(); return 1;
        }
    }

    if (evalMode) {
        if (gtFile.empty() || bedFile.empty()) {
            std::cerr << "Error: --evaluate requires -g and -b\n";
            printUsage();
            return 1;
        }
        evaluate(gtFile, bedFile);
    } else {
        generate(params);
    }

    return 0;
}
