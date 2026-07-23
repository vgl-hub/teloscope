#include <iostream>
#include <algorithm>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <tuple>
#include <fstream>
#include <stdexcept> // jack: std::runtime_error
#include <cctype>
#include <cstring>

#include "log.h"
#include "global.h"
#include "uid-generator.h"
#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "stream-obj.h"
#include "fastx.h"
#include "gfa-lines.h"
#include "gfa.h"
#include "input-gfa.h"
#include "threadpool.h"
#include "teloscope.h"
#include "output.h"
#include "input.h"
#include "read-filter.h"

namespace {

struct PendingTelomereAnnotation {
    std::string dedupeKey;
    std::string header;
    bool atStart = false;
    uint32_t blockLen = 0;
    char segOrient = '+';
};

struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string plus;
    std::string quality;
};

struct FastqChunkResult {
    std::string output;
    uint64_t scanned = 0;
    uint64_t passed = 0;
};

std::vector<Tag> makeSyntheticTelomereTags(uint32_t blockLen) {
    return {
        Tag{'i', "LN", "6"},
        Tag{'i', "RC", "6000"},
        Tag{'i', "TL", std::to_string(blockLen)}
    };
}

void queueAnnotation(std::vector<PendingTelomereAnnotation> &annotations,
                     const std::string &dedupeKey,
                     const std::string &header,
                     bool atStart,
                     uint32_t blockLen,
                     char segOrient) {
    for (auto &ann : annotations) {
        if (ann.dedupeKey == dedupeKey) {
            ann.blockLen = std::max(ann.blockLen, blockLen);
            return;
        }
    }

    annotations.push_back({dedupeKey, header, atStart, blockLen, segOrient});
}

void appendTelomereConnection(InSequences &inSequences,
                              unsigned int telomereUid,
                              unsigned int segmentUid,
                              char segmentOrient) {
    // Telomere cap is a direct adjacency: an L link with 0M overlap, not a J gap.
    InEdge edge;
    edge.newEdge(0, telomereUid, segmentUid, '+', segmentOrient, "0M", "",
                 std::vector<Tag>{Tag{'i', "RC", "0"}});
    inSequences.appendEdge(edge);
}

[[noreturn]] void fastqExitFailure() {
    threadPool.join();
    exit(EXIT_FAILURE);
}

size_t logicalLineLength(const std::string &line) {
    if (!line.empty() && line.back() == '\r') {
        return line.size() - 1;
    }
    return line.size();
}

[[noreturn]] void fastqInputError(uint64_t recordNumber, const std::string &message) {
    if (recordNumber > 0) {
        fprintf(stderr, "Error: FASTQ record %" PRIu64 ": %s.\n",
                recordNumber, message.c_str());
    } else {
        fprintf(stderr, "Error: %s.\n", message.c_str());
    }
    fastqExitFailure();
}

bool readFastqRecord(std::istream &stream, FastqRecord &record, uint64_t recordNumber) {
    // Skip blank lines (including a lone '\r') before a header so a trailing newline
    // or a stray blank line does not abort an otherwise valid stream.
    do {
        if (!std::getline(stream, record.header)) {
            return false;
        }
    } while (logicalLineLength(record.header) == 0);
    if (!std::getline(stream, record.sequence) ||
        !std::getline(stream, record.plus) ||
        !std::getline(stream, record.quality)) {
        fastqInputError(recordNumber, "truncated FASTQ record");
    }

    if (record.header.empty() || record.header.front() != '@') {
        fastqInputError(recordNumber, "expected header line starting with '@'");
    }
    if (record.plus.empty() || record.plus.front() != '+') {
        fastqInputError(recordNumber, "expected separator line starting with '+'");
    }
    if (logicalLineLength(record.sequence) != logicalLineLength(record.quality)) {
        fastqInputError(recordNumber, "sequence and quality length differ");
    }

    return true;
}

void appendFastqRecord(std::string &out, const FastqRecord &record) {
    out += record.header;
    out += '\n';
    out += record.sequence;
    out += '\n';
    out += record.plus;
    out += '\n';
    out += record.quality;
    out += '\n';
}

[[noreturn]] void sequenceFilterError(const std::string &message) {
    fprintf(stderr, "Error: %s\n", message.c_str());
    threadPool.join();
    exit(EXIT_FAILURE);
}

std::string trimFilterLine(const std::string &value) {
    const size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string sequenceFilterId(const std::string &header) {
    const size_t firstWhitespace = header.find_first_of(" \t\r\n\f\v");
    return header.substr(0, firstWhitespace);
}

bool hasCaseInsensitiveSuffix(const std::string &value, const std::string &suffix) {
    if (value.size() < suffix.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), value.rbegin(),
        [](unsigned char left, unsigned char right) {
            return std::tolower(left) == std::tolower(right);
        });
}

bool isGfaAssemblyPath(const std::string &path) {
    return hasCaseInsensitiveSuffix(path, ".gfa") ||
           hasCaseInsensitiveSuffix(path, ".gfa.gz") ||
           hasCaseInsensitiveSuffix(path, ".gfa2") ||
           hasCaseInsensitiveSuffix(path, ".gfa2.gz");
}

bool readGzipLine(gzFile input, std::vector<char> &buffer, std::string &line,
                  const std::string &path) {
    line.clear();
    while (true) {
        char *chunk = gzgets(input, buffer.data(), static_cast<int>(buffer.size()));
        if (chunk == nullptr) {
            if (gzeof(input)) return !line.empty();
            int errorNumber = Z_OK;
            const char *errorMessage = gzerror(input, &errorNumber);
            sequenceFilterError("Could not read assembly input '" + path + "': " +
                                (errorMessage ? errorMessage : "zlib error") + ".");
        }
        const size_t chunkLength = std::strlen(chunk);
        line.append(chunk, chunkLength);
        if (chunkLength > 0 && line.back() == '\n') {
            line.pop_back();
            return true;
        }
        if (gzeof(input)) return true;
    }
}

void validateFilteredGfaInput(const UserInputTeloscope &input) {
    if (hasCaseInsensitiveSuffix(input.inSequence, ".gfa2") ||
        hasCaseInsensitiveSuffix(input.inSequence, ".gfa2.gz")) {
        sequenceFilterError("Assembly record filters do not support GFA2; use GFA1 P paths or a pathless GFA1 graph.");
    }

    gzFile stream = gzopen(input.inSequence.c_str(), "rb");
    if (stream == nullptr) {
        sequenceFilterError("Could not open assembly input '" + input.inSequence + "'.");
    }
    gzbuffer(stream, 1U << 20);
    std::vector<char> buffer(1U << 16);
    std::string validationError;
    uint64_t lineNumber = 0;
    for (std::string line; readGzipLine(stream, buffer, line, input.inSequence); ) {
        lineNumber++;
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty() || line.front() == '#') continue;
        if (line.rfind("H\t", 0) == 0 && line.find("\tVN:Z:2") != std::string::npos) {
            validationError = "Assembly record filters do not support GFA2 at line " +
                              std::to_string(lineNumber) +
                              "; use GFA1 P paths or a pathless GFA1 graph.";
            break;
        }
        if (line.size() < 2 || line[1] != '\t') {
            validationError = "Assembly record filters found a malformed or unsupported GFA record at line " +
                              std::to_string(lineNumber) + ".";
            break;
        }
        const char recordType = line.front();
        if (recordType == 'O' || recordType == 'U' || recordType == 'E' ||
            recordType == 'G' || recordType == 'F') {
            validationError = std::string("Assembly record filters do not support GFA2 record type '") +
                              recordType + "' at line " + std::to_string(lineNumber) +
                              "; use GFA1 P paths or a pathless GFA1 graph.";
            break;
        }
        if (recordType == 'W') {
            validationError = "Assembly record filters do not support GFA1 W walks at line " +
                              std::to_string(lineNumber) +
                              "; use GFA1 P paths or a pathless GFA1 graph.";
            break;
        }
        if (recordType == 'C') {
            validationError = "Assembly record filters do not support GFA1 C containment records at line " +
                              std::to_string(lineNumber) + ".";
            break;
        }
        if (recordType == 'S') {
            const size_t secondTab = line.find('\t', 2);
            const size_t thirdTab = secondTab == std::string::npos
                ? std::string::npos : line.find('\t', secondTab + 1);
            if (thirdTab != std::string::npos) {
                const std::string lengthField = line.substr(secondTab + 1,
                                                            thirdTab - secondTab - 1);
                if (!lengthField.empty() &&
                    std::all_of(lengthField.begin(), lengthField.end(),
                                [](unsigned char c) { return std::isdigit(c); })) {
                    validationError = "Assembly record filters do not support GFA2 segment records at line " +
                                      std::to_string(lineNumber) +
                                      "; use GFA1 P paths or a pathless GFA1 graph.";
                    break;
                }
            }
        }
        if (recordType != 'H' && recordType != 'S' && recordType != 'L' &&
            recordType != 'J' && recordType != 'P') {
            validationError = std::string("Assembly record filters do not support GFA record type '") +
                              recordType + "' at line " + std::to_string(lineNumber) + ".";
            break;
        }
    }
    const int closeResult = gzclose(stream);
    if (closeResult != Z_OK) {
        sequenceFilterError("Could not close assembly input '" + input.inSequence + "'.");
    }
    if (!validationError.empty()) sequenceFilterError(validationError);
}

void loadNormalizedFastaAssembly(UserInputTeloscope &input, InSequences &sequences) {
    gzFile gzipInput = nullptr;
    std::istream *plainInput = nullptr;
    if (input.inSequence.empty()) {
        plainInput = &std::cin;
    } else {
        gzipInput = gzopen(input.inSequence.c_str(), "rb");
        if (gzipInput == nullptr) {
            sequenceFilterError("Could not open assembly input '" + input.inSequence + "'.");
        }
        gzbuffer(gzipInput, 1U << 20);
    }

    std::vector<char> gzipBuffer(1U << 16);
    auto readLine = [&](std::string &line) {
        if (plainInput != nullptr) return static_cast<bool>(std::getline(*plainInput, line));
        return readGzipLine(gzipInput, gzipBuffer, line, input.inSequence);
    };

    uint32_t sequencePosition = 0;
    std::unordered_set<std::string> seenIds;
    std::string primaryId;
    std::string comment;
    std::string *sequence = nullptr;
    auto appendRecord = [&]() {
        if (sequence == nullptr) return;
        if (sequence->empty()) {
            delete sequence;
            sequenceFilterError("FASTA record '" + primaryId + "' has no sequence.");
        }
        Sequence *record = new Sequence{primaryId, comment, sequence, nullptr, sequencePosition++};
        sequences.appendSequence(record, input.hc_cutoff);
        sequence = nullptr;
    };

    bool firstLine = true;
    for (std::string line; readLine(line); ) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (firstLine && line.size() >= 3 &&
            static_cast<unsigned char>(line[0]) == 0xef &&
            static_cast<unsigned char>(line[1]) == 0xbb &&
            static_cast<unsigned char>(line[2]) == 0xbf) {
            line.erase(0, 3);
        }
        firstLine = false;

        if (!line.empty() && line.front() == '>') {
            appendRecord();
            const std::string headerLine = line.substr(1);
            primaryId = sequenceFilterId(headerLine);
            if (primaryId.empty()) {
                sequenceFilterError("FASTA input contains an empty primary sequence ID.");
            }
            if (!seenIds.insert(primaryId).second) {
                sequenceFilterError("Input contains duplicate primary sequence ID: '" + primaryId + "'.");
            }
            const size_t firstWhitespace = headerLine.find_first_of(" \t\r\n\f\v");
            comment = firstWhitespace == std::string::npos
                ? "" : trimFilterLine(headerLine.substr(firstWhitespace + 1));
            sequence = new std::string;
        } else {
            if (sequence == nullptr) {
                sequenceFilterError("Assembly record filters require FASTA input or a recognized GFA file.");
            }
            sequence->append(line);
        }
    }

    appendRecord();
    if (sequencePosition == 0) sequenceFilterError("Assembly input is empty.");
    if (gzipInput != nullptr && gzclose(gzipInput) != Z_OK) {
        sequenceFilterError("Could not close assembly input '" + input.inSequence + "'.");
    }

    jobWait(threadPool);
    sequences.updateStats();
}

bool parseUnsignedCoordinate(const std::string &value, uint64_t &coordinate) {
    if (value.empty() ||
        !std::all_of(value.begin(), value.end(), [](unsigned char c) { return std::isdigit(c); })) {
        return false;
    }
    try {
        size_t parsed = 0;
        const unsigned long long result = std::stoull(value, &parsed);
        if (parsed != value.size()) return false;
        coordinate = static_cast<uint64_t>(result);
        return true;
    } catch (...) {
        return false;
    }
}

struct SequenceSelection {
    std::unordered_set<std::string> names;
    uint64_t selectedCount = 0;
};

class SequenceSelector {
    std::unordered_set<std::string> includeIds;
    std::unordered_set<std::string> excludeIds;
    std::vector<std::string> includePrefixes;
    std::vector<std::string> excludePrefixes;
    bool active = false;

    static bool startsWith(const std::string &value, const std::string &prefix) {
        return value.size() >= prefix.size() &&
               value.compare(0, prefix.size(), prefix) == 0;
    }

    static void deduplicatePrefixes(std::vector<std::string> &prefixes) {
        std::unordered_set<std::string> seen;
        std::vector<std::string> unique;
        unique.reserve(prefixes.size());
        for (const std::string &prefix : prefixes) {
            if (seen.insert(prefix).second) unique.push_back(prefix);
        }
        prefixes.swap(unique);
    }

    static void loadSelectorFiles(const std::vector<std::string> &files,
                                  std::unordered_set<std::string> &ids,
                                  const char *optionName) {
        for (const std::string &path : files) {
            std::ifstream input(path);
            if (!input.is_open()) {
                sequenceFilterError(std::string("Could not open ") + optionName + " file '" + path + "'.");
            }

            uint64_t idLines = 0;
            uint64_t lineNumber = 0;
            for (std::string raw; std::getline(input, raw); ) {
                lineNumber++;
                if (lineNumber == 1 && raw.size() >= 3 &&
                    static_cast<unsigned char>(raw[0]) == 0xef &&
                    static_cast<unsigned char>(raw[1]) == 0xbb &&
                    static_cast<unsigned char>(raw[2]) == 0xbf) {
                    raw.erase(0, 3);
                }
                const std::string line = trimFilterLine(raw);
                if (line.empty() || line.front() == '#') continue;

                std::istringstream fieldsStream(line);
                std::vector<std::string> fields;
                for (std::string field; fieldsStream >> field; ) fields.push_back(field);
                if (fields.empty()) continue;
                if (fields[0] == "track" || fields[0] == "browser") continue;
                if (fields.size() == 2) {
                    sequenceFilterError(path + ":" + std::to_string(lineNumber) +
                                        " must contain either one ID column or at least three BED columns.");
                }
                if (fields.size() >= 3) {
                    uint64_t begin = 0, end = 0;
                    if (!parseUnsignedCoordinate(fields[1], begin) ||
                        !parseUnsignedCoordinate(fields[2], end) || begin > end) {
                        sequenceFilterError(path + ":" + std::to_string(lineNumber) +
                                            " has invalid BED start/end coordinates.");
                    }
                }

                ids.insert(fields[0]);
                idLines++;
            }

            if (idLines == 0) {
                sequenceFilterError(std::string(optionName) + " file '" + path +
                                    "' contains no sequence IDs.");
            }
        }
    }

    static bool matchesAnyPrefix(const std::string &name,
                                 const std::vector<std::string> &prefixes) {
        return std::any_of(prefixes.begin(), prefixes.end(),
                           [&](const std::string &prefix) { return startsWith(name, prefix); });
    }

    static std::string describeUnmatched(const std::vector<std::string> &values) {
        constexpr size_t displayLimit = 10;
        std::ostringstream message;
        const size_t shown = std::min(values.size(), displayLimit);
        for (size_t i = 0; i < shown; ++i) {
            if (i > 0) message << ", ";
            message << "'" << values[i] << "'";
        }
        if (values.size() > shown) message << " (and " << (values.size() - shown) << " more)";
        return message.str();
    }

    void validateSelectors(const std::vector<std::string> &names,
                           const std::string &domainLabel) const {
        const std::unordered_set<std::string> available(names.begin(), names.end());
        std::vector<std::string> unmatchedIds;
        std::vector<std::string> unmatchedPrefixes;

        for (const std::string &id : includeIds)
            if (available.count(id) == 0) unmatchedIds.push_back(id);
        for (const std::string &id : excludeIds)
            if (available.count(id) == 0) unmatchedIds.push_back(id);

        auto validatePrefixList = [&](const std::vector<std::string> &prefixes) {
            for (const std::string &prefix : prefixes) {
                const bool matched = std::any_of(names.begin(), names.end(),
                    [&](const std::string &name) { return startsWith(name, prefix); });
                if (!matched) unmatchedPrefixes.push_back(prefix);
            }
        };
        validatePrefixList(includePrefixes);
        validatePrefixList(excludePrefixes);

        if (!unmatchedIds.empty()) {
            std::sort(unmatchedIds.begin(), unmatchedIds.end());
            unmatchedIds.erase(std::unique(unmatchedIds.begin(), unmatchedIds.end()),
                               unmatchedIds.end());
            sequenceFilterError("Sequence filter ID(s) matched no input " + domainLabel + ": " +
                                describeUnmatched(unmatchedIds) + ".");
        }
        if (!unmatchedPrefixes.empty()) {
            std::sort(unmatchedPrefixes.begin(), unmatchedPrefixes.end());
            unmatchedPrefixes.erase(std::unique(unmatchedPrefixes.begin(), unmatchedPrefixes.end()),
                                    unmatchedPrefixes.end());
            sequenceFilterError("Sequence filter prefix(es) matched no input " + domainLabel + ": " +
                                describeUnmatched(unmatchedPrefixes) + ".");
        }
    }

public:
    explicit SequenceSelector(const UserInputTeloscope &input)
        : includePrefixes(input.includePrefixes), excludePrefixes(input.excludePrefixes),
          active(input.sequenceFilterActive) {
        loadSelectorFiles(input.includeBedFiles, includeIds, "--include-bed");
        loadSelectorFiles(input.excludeBedFiles, excludeIds, "--exclude-bed");
        deduplicatePrefixes(includePrefixes);
        deduplicatePrefixes(excludePrefixes);
    }

    SequenceSelection select(const std::vector<std::string> &candidateNames,
                             const std::string &domainLabel) const {
        SequenceSelection selection;
        selection.names.reserve(candidateNames.size());

        if (active) {
            std::unordered_set<std::string> uniqueNames;
            std::vector<std::string> duplicateNames;
            for (const std::string &name : candidateNames) {
                if (name.empty()) {
                    sequenceFilterError("Input contains an empty primary sequence ID.");
                }
                if (!uniqueNames.insert(name).second) duplicateNames.push_back(name);
            }
            if (!duplicateNames.empty()) {
                std::sort(duplicateNames.begin(), duplicateNames.end());
                duplicateNames.erase(std::unique(duplicateNames.begin(), duplicateNames.end()),
                                     duplicateNames.end());
                sequenceFilterError("Input contains duplicate primary sequence ID(s): " +
                                    describeUnmatched(duplicateNames) + ".");
            }
            validateSelectors(candidateNames, domainLabel);
        }

        const bool hasIncludes = !includeIds.empty() || !includePrefixes.empty();
        for (const std::string &name : candidateNames) {
            const bool included = !hasIncludes || includeIds.count(name) != 0 ||
                                  matchesAnyPrefix(name, includePrefixes);
            const bool excluded = excludeIds.count(name) != 0 ||
                                  matchesAnyPrefix(name, excludePrefixes);
            if (included && !excluded) {
                selection.names.insert(name);
                selection.selectedCount++;
            }
        }

        if (active && selection.selectedCount == 0) {
            sequenceFilterError("Sequence filters excluded all input " + domainLabel + ".");
        }
        return selection;
    }
};

} // namespace


void Input::load(UserInputTeloscope userInput) {
    
    this->userInput = userInput;
    
}


void Input::read(InSequences &inSequences) {
    SequenceSelector selector(userInput);
    const bool isGfa = isGfaAssemblyPath(userInput.inSequence);
    if (isGfa && userInput.sequenceFilterActive) validateFilteredGfaInput(userInput);
    if (isGfa || !userInput.sequenceFilterActive) {
        loadGenome(userInput, inSequences);
    } else {
        loadNormalizedFastaAssembly(userInput, inSequences);
    }
    lg.verbose("Finished loading genome assembly");

    std::vector<InSegment*> *inSegments = inSequences.getInSegments();
    std::vector<InGap> *inGaps = inSequences.getInGaps();

    std::vector<InPath> inPaths = inSequences.getInPaths();
    if (!isGfa && userInput.sequenceFilterActive) {
        // Normalize FASTA primary IDs at the first ASCII whitespace character.
        for (InPath &path : inPaths) path.setHeader(sequenceFilterId(path.getHeader()));
    }

    const bool filterSegments = isGfa && inPaths.empty();
    std::vector<std::string> candidateNames;
    if (filterSegments) {
        candidateNames.reserve(inSegments->size());
        for (InSegment *segment : *inSegments)
            candidateNames.push_back(sequenceFilterId(segment->getSeqHeader()));
    } else {
        candidateNames.reserve(inPaths.size());
        for (InPath &path : inPaths)
            candidateNames.push_back(sequenceFilterId(path.getHeader()));
    }

    const std::string domainLabel = filterSegments ? "segments" : "paths";
    const SequenceSelection selection = selector.select(candidateNames, domainLabel);
    if (userInput.sequenceFilterActive) {
        userInput.filterInputCount = candidateNames.size();
        userInput.filterSelectedCount = selection.selectedCount;
        fprintf(stderr, "Sequence filter: selected %" PRIu64 " of %" PRIu64 " %s.\n",
                userInput.filterSelectedCount, userInput.filterInputCount, domainLabel.c_str());
    }
    if (!filterSegments) {
        inPaths.erase(std::remove_if(inPaths.begin(), inPaths.end(),
            [&](InPath &path) {
                return selection.names.count(sequenceFilterId(path.getHeader())) == 0;
            }), inPaths.end());
    }

    Teloscope teloscope(userInput);

    // GFA-based annotation
    if (isGfa) {

        // One scan job per terminal end. Resolve every job before queuing any:
        // worker threads append telomere nodes to the live segment vector, so
        // reading or iterating it while jobs run would race that growth.
        struct TeloJob { InSegment* seg; char orient; bool isFirst; bool pathAware; };
        std::vector<TeloJob> jobs;

        if (!inPaths.empty()) {
            // Unique (segment, orientation, isFirst) terminal ends: one job per node.
            std::set<std::tuple<unsigned int, char, bool>> terminalEnds;
            for (InPath& path : inPaths) {
                std::vector<PathComponent> components = path.getComponents();

                // First SEGMENT component => contig start (isFirst = true)
                for (const auto& comp : components) {
                    if (comp.componentType == SEGMENT && comp.orientation != '0') {
                        terminalEnds.emplace(comp.id, comp.orientation, true);
                        break;
                    }
                }

                // Last SEGMENT component => contig end (isFirst = false)
                for (auto it = components.rbegin(); it != components.rend(); ++it) {
                    if (it->componentType == SEGMENT && it->orientation != '0') {
                        terminalEnds.emplace(it->id, it->orientation, false);
                        break;
                    }
                }
            }

            for (const auto& end : terminalEnds)
                jobs.push_back({inSequences.getInSegment(std::get<0>(end)),
                                std::get<1>(end), std::get<2>(end), true});
        } else {
            // Fallback: path-less GFA, scan all segments (implicit + orientation)
            for (InSegment* inSegment : *inSegments) {
                if (selection.names.count(sequenceFilterId(inSegment->getSeqHeader())) != 0)
                    jobs.push_back({inSegment, '+', false, false});
            }
        }

        // Count missing sequence now, while the segment vector is still stable.
        size_t segmentsScanned = jobs.size(), segmentsNoSeq = 0;
        for (const TeloJob& job : jobs)
            if (job.seg->getInSequencePtr() == NULL) ++segmentsNoSeq;

        // No DNA means nothing to scan: say so instead of silently annotating nothing.
        if (segmentsNoSeq > 0)
            fprintf(stderr, "Warning: %zu of %zu GFA segment(s) had no sequence (*); skipped for telomere annotation.\n",
                    segmentsNoSeq, segmentsScanned);

        for (const TeloJob& job : jobs) {
            InSegment* seg = job.seg;
            char orient = job.orient;
            bool isFirst = job.isFirst;
            if (job.pathAware)
                threadPool.queueJob([seg, &inSequences, &teloscope, orient, isFirst]() {
                    return teloscope.walkSegmentForPath(seg, inSequences, orient, isFirst);
                });
            else
                threadPool.queueJob([seg, &inSequences, &teloscope]() {
                    return teloscope.walkSegment(seg, inSequences);
                });
        }

        lg.verbose("Waiting for telomere annotation jobs to complete");
        jobWait(threadPool);
        lg.verbose("\nAll telomere annotation jobs completed.");

        // Write annotated GFA
        Report report;
        std::string outGfa = userInput.outRoute + "/" + userInput.inSequenceName + ".telo.annotated.gfa";
        report.writeToStream(inSequences, outGfa, userInput);
        lg.verbose("Annotated GFA written to " + outGfa);

        // Companion colours CSV: every telomere node green (#008000), read from the graph.
        std::string outColors = userInput.outRoute + "/" + userInput.inSequenceName + ".telo.annotated.colors.csv";
        std::ofstream colorsStream(outColors);
        if (colorsStream) {
            colorsStream << "node\tcolor\n";
            for (InSegment* seg : *inSequences.getInSegments()) {
                const std::string& name = seg->getSeqHeader();
                if (name.rfind("telomere_", 0) == 0)
                    colorsStream << name << "\t#008000\n";
            }
            lg.verbose("Telomere colours CSV written to " + outColors);
        } else {
            fprintf(stderr, "Warning: could not write telomere colours CSV to %s.\n", outColors.c_str());
        }
        return;
    }

    // path-based annotation
    for (InPath& inPath : inPaths) {
        InPath* pathPtr = &inPath;
        threadPool.queueJob([pathPtr, inSegments, inGaps, &teloscope]() {
            return teloscope.walkPath(pathPtr, *inSegments, *inGaps);
        });
    }
    lg.verbose("Waiting for jobs to complete");
    jobWait(threadPool);
    lg.verbose("\nAll jobs completed.");

    teloscope.sortBySeqPos();
    lg.verbose("\nPaths sorted by original position.");

    teloscope.handleBEDFile();
    lg.verbose("\nReport and BED/BEDgraph files generated.");
}


void Input::readFastqSubset(std::ostream &out) {
    StreamObj streamObj;
    std::shared_ptr<std::istream> stream = streamObj.openStream(userInput, 'f');
    if (!stream) {
        fprintf(stderr, "Error: Stream not successful: %s.\n", userInput.inSequence.c_str());
        fastqExitFailure();
    }

    int first = stream->peek();
    if (first == EOF) {
        fastqInputError(0, "FASTQ input is empty");
    }
    if (first != '@') {
        fastqInputError(0, "FASTQ input must start with '@'");
    }

    const uint32_t threads = std::max<uint32_t>(1, threadPool.totalThreads());
    const size_t recordsPerBatch = std::min<size_t>(
        2048, std::max<size_t>(256, static_cast<size_t>(threads) * 32));

    std::vector<FastqRecord> batch;
    batch.reserve(recordsPerBatch);

    uint64_t recordNumber = 0;
    uint64_t totalReads = 0;
    uint64_t passedReads = 0;

    auto processBatch = [&]() {
        if (batch.empty()) {
            return;
        }

        const size_t chunkCount = std::min<size_t>(threads, batch.size());
        const size_t chunkSize = (batch.size() + chunkCount - 1) / chunkCount;
        std::vector<FastqChunkResult> results(chunkCount);

        for (size_t chunk = 0; chunk < chunkCount; ++chunk) {
            const size_t start = chunk * chunkSize;
            const size_t end = std::min(batch.size(), start + chunkSize);
            if (start >= end) {
                continue;
            }

            threadPool.queueJob([&, chunk, start, end]() {
                ReadTelomereFilter filter(userInput);
                FastqChunkResult &result = results[chunk];
                result.scanned = end - start;

                for (size_t i = start; i < end; ++i) {
                    if (filter.matches(batch[i].sequence)) {
                        appendFastqRecord(result.output, batch[i]);
                        result.passed++;
                    }
                }

                return true;
            });
        }

        jobWait(threadPool);

        for (const auto &result : results) {
            if (!result.output.empty()) {
                out.write(result.output.data(), static_cast<std::streamsize>(result.output.size()));
            }
            totalReads += result.scanned;
            passedReads += result.passed;
        }

        if (!out.good()) {
            fprintf(stderr, "Error: failed while writing FASTQ subset to stdout.\n");
            fastqExitFailure();
        }

        batch.clear();
    };

    while (true) {
        FastqRecord record;
        if (!readFastqRecord(*stream, record, recordNumber + 1)) {
            break;
        }
        recordNumber++;
        batch.push_back(std::move(record));

        if (batch.size() == recordsPerBatch) {
            processBatch();
        }
    }

    processBatch();
    out.flush();

    fprintf(stderr, "FASTQ subset: kept %" PRIu64 " of %" PRIu64 " reads.\n",
            passedReads, totalReads);
}


bool Teloscope::walkSegment(InSegment* segment, InSequences& inSequences) {
    Log threadLog;
    threadLog.add("\n\tWalking segment:\t" + segment->getSeqHeader());

    std::string sequence = segment->getInSequence(0, 0);
    unmaskSequence(sequence);

    SegmentData segmentData = scanSegment(sequence, 0, true); // tipsOnly = true for GFA segments

    std::vector<PendingTelomereAnnotation> annotations;

    for (const TelomereBlock& block : segmentData.terminalBlocks) {
        uint64_t distToStart = block.start;
        uint64_t distToEnd   = sequence.size() - (block.start + block.blockLen);
        bool atStart = distToStart <= distToEnd;

        // edge orientation: + = start, - = end
        char segOrient = atStart ? '+' : '-';
        // Path-less segments are scanned as '+', so encode '+' in the node name.
        const std::string header = "telomere_"
                        + segment->getSeqHeader()
                        + "+"
                        + (atStart ? "_start" : "_end");
        queueAnnotation(annotations, header, header, atStart, block.blockLen, segOrient);
    }

    for (auto& ann : annotations) {
        Sequence teloSeq{ann.header, "", new std::string("*")};
        inSequences.traverseInSegmentWrapper(&teloSeq, makeSyntheticTelomereTags(ann.blockLen));

        unsigned int teloUid;
        {
            std::lock_guard<std::mutex> lck(mtx);
            teloUid = inSequences.getHash1()->at(ann.header);
        }

        appendTelomereConnection(inSequences, teloUid, segment->getuId(), ann.segOrient);
    }

    threadLog.add("\tCompleted walking segment:\t" + segment->getSeqHeader());
    {
        std::lock_guard<std::mutex> lck(mtx);
        logs.push_back(threadLog);
    }

    return true;
}


bool Teloscope::walkSegmentForPath(InSegment* segment, InSequences& inSequences,
                                   char pathOrient, bool isFirst) {
    Log threadLog;
    threadLog.add("\n\tWalking segment (path-aware):\t" + segment->getSeqHeader());

    std::string sequence = segment->getInSequence(0, 0);
    unmaskSequence(sequence);

    SegmentData segmentData = scanSegment(sequence, 0, true);

    // Physical end holding the path-terminal tip (start when isFirst == (orient=='+')).
    bool scanStart = (isFirst == (pathOrient == '+'));

    std::vector<PendingTelomereAnnotation> annotations;

    for (const TelomereBlock& block : segmentData.terminalBlocks) {
        uint64_t distToStart = block.start;
        uint64_t distToEnd   = sequence.size() - (block.start + block.blockLen);
        bool blockAtStart = distToStart <= distToEnd;

        // Only keep the block at the path-terminal end
        if (blockAtStart != scanStart) continue;

        std::string posLabel = isFirst ? "start" : "end";
        std::string header = "telomere_"
                        + segment->getSeqHeader()
                        + std::string(1, pathOrient)
                        + "_" + posLabel;

        // Edge orientation rule: START = pathOrient, END = flip(pathOrient)
        char edgeOrient = isFirst ? pathOrient : (pathOrient == '+' ? '-' : '+');
        queueAnnotation(annotations, header, header, blockAtStart, block.blockLen, edgeOrient);
    }

    for (auto& ann : annotations) {
        // Terminal ends were de-duplicated upstream, so each header is created once.
        Sequence teloSeq{ann.header, "", new std::string("*")};
        inSequences.traverseInSegmentWrapper(&teloSeq, makeSyntheticTelomereTags(ann.blockLen));

        unsigned int teloUid;
        {
            std::lock_guard<std::mutex> lck(mtx);
            teloUid = inSequences.getHash1()->at(ann.header);
        }

        appendTelomereConnection(inSequences, teloUid, segment->getuId(), ann.segOrient);
    }

    threadLog.add("\tCompleted walking segment (path-aware):\t" + segment->getSeqHeader());
    {
        std::lock_guard<std::mutex> lck(mtx);
        logs.push_back(threadLog);
    }

    return true;
}


bool Teloscope::walkPath(InPath* path, std::vector<InSegment*> &inSegments, std::vector<InGap> &inGaps) {
    Log threadLog;
    uint64_t absPos = 0;
    unsigned int cUId = 0, gapLen = 0, seqPos = path->getSeqPos();
    std::vector<PathComponent> pathComponents = path->getComponents();

    threadLog.add("\n\tWalking path:\t" + path->getHeader());
    std::string header = path->getHeader();
    eraseChar(header, '\r');

    // Initialize PathData for this path
    PathData pathData;
    pathData.seqPos = seqPos;
    pathData.header = header;
    pathData.pathSize = path->getLen();
    // pathData.windows.reserve(inSegments.size()); NumWindows = ceil((L - W) / S) + 1

    // index for O(1) lookups
    std::unordered_map<unsigned int, InSegment*> segmentIndex;
    segmentIndex.reserve(inSegments.size());
    for (auto* seg : inSegments) segmentIndex[seg->getuId()] = seg;

    std::unordered_map<unsigned int, InGap*> gapIndex;
    gapIndex.reserve(inGaps.size());
    for (auto& gap : inGaps) gapIndex[gap.getuId()] = &gap;

    for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
        cUId = component->id;

        if (component->componentType == SEGMENT) {
            auto inSegment = segmentIndex.find(cUId)->second;
            std::string sequence = inSegment->getInSequence(component->start, component->end);
            unmaskSequence(sequence);
            
            if (component->orientation == '+') {
                SegmentData segmentData = scanSegment(sequence, absPos, userInput.ultraFastMode);

                // Collect window data
                pathData.windows.insert(
                    pathData.windows.end(),
                    std::make_move_iterator(segmentData.windows.begin()),
                    std::make_move_iterator(segmentData.windows.end())
                );

                // Collect blocks
                pathData.terminalBlocks.insert(
                    pathData.terminalBlocks.end(),
                    std::make_move_iterator(segmentData.terminalBlocks.begin()),
                    std::make_move_iterator(segmentData.terminalBlocks.end())
                );

                pathData.interstitialBlocks.insert(
                    pathData.interstitialBlocks.end(),
                    std::make_move_iterator(segmentData.interstitialBlocks.begin()),
                    std::make_move_iterator(segmentData.interstitialBlocks.end())
                );

                // Collect matches
                pathData.canonicalMatches.insert(
                    pathData.canonicalMatches.end(),
                    std::make_move_iterator(segmentData.canonicalMatches.begin()),
                    std::make_move_iterator(segmentData.canonicalMatches.end())
                );

                pathData.nonCanonicalMatches.insert(
                    pathData.nonCanonicalMatches.end(),
                    std::make_move_iterator(segmentData.nonCanonicalMatches.begin()),
                    std::make_move_iterator(segmentData.nonCanonicalMatches.end())
                );

            } else {
            }
            
            absPos += sequence.size();
            
        }else if (component->componentType == GAP){
            
            auto inGap = gapIndex.find(cUId)->second;
            gapLen = inGap->getDist(component->start - component->end);

            pathData.gapInfos.push_back({absPos, static_cast<uint32_t>(gapLen)});
            absPos += gapLen;
            
        } else {
        } // need to handle edges, cigars etc
        
    }

    // Filter blocks
    labelTerminalBlocks(pathData.terminalBlocks, static_cast<uint16_t>(pathData.gapInfos.size()),
                        pathData.terminalLabel, pathData.scaffoldType,
                        pathData.pathSize, userInput.terminalLimit);
    threadLog.add("\tCompleted walking path:\t" + path->getHeader());

    std::lock_guard<std::mutex> lck(mtx);
    allPathData.push_back(std::move(pathData));
    logs.push_back(threadLog);

    return true;
}
