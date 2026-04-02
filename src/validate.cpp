/*
USAGE:
test <path to test folder or files>

EXAMPLE:
build/bin/teloscope-validate validateFiles
build/bin/teloscope-validate validateFiles/random1.fasta0.tst
*/

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits.h>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unistd.h>
#include <vector>

#ifndef _WIN32
#include <sys/wait.h>
#endif

#include <validate.h>

namespace fs = std::filesystem;

bool printCommand = false;
bool pass = true;

void printFAIL(const char *m1="", const char *m2="", const char *m3="", const char *m4="") {
    pass = false;
    std::cout << "\033[0;31mFAIL\033[0m " << m1 << " " << m2 << " " << m3 << " " << m4 << std::endl;
}

void printPASS(const char *m1="", const char *m2="", const char *m3="", const char *m4="") {
    std::cout << "\033[0;32mPASS\033[0m " << m1 << " " << m2 << " " << m3 << " " << m4 << std::endl;
}

namespace {

std::string trim(const std::string &line) {
    const auto begin = line.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos)
        return "";
    const auto end = line.find_last_not_of(" \t\r\n");
    return line.substr(begin, end - begin + 1);
}

bool startsWith(const std::string &text, const std::string &prefix) {
    return text.rfind(prefix, 0) == 0;
}

std::vector<std::string> split(const std::string &text, char delim) {
    std::vector<std::string> fields;
    std::stringstream stream(text);
    std::string field;
    while (std::getline(stream, field, delim))
        fields.push_back(field);
    return fields;
}

std::string joinTail(const std::vector<std::string> &parts, size_t start, char delim) {
    std::ostringstream joined;
    for (size_t i = start; i < parts.size(); ++i) {
        if (i != start)
            joined << delim;
        joined << parts[i];
    }
    return joined.str();
}

std::string stripCarriageReturn(std::string line) {
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    return line;
}

std::string readFile(const fs::path &path) {
    std::ifstream stream(path);
    std::ostringstream content;
    content << stream.rdbuf();
    return content.str();
}

int normalizeExitCode(int status) {
    if (status == -1)
        return -1;
#ifdef _WIN32
    return status;
#else
    if (WIFEXITED(status))
        return WEXITSTATUS(status);
    if (WIFSIGNALED(status))
        return 128 + WTERMSIG(status);
    return status;
#endif
}

std::string shellQuote(const std::string &text) {
    return "\"" + text + "\"";
}

std::string sanitizeName(const std::string &text) {
    std::string out = text;
    for (char &c : out) {
        if (!std::isalnum(static_cast<unsigned char>(c)))
            c = '_';
    }
    return out;
}

std::string replaceAll(std::string text, const std::string &needle, const std::string &replacement) {
    size_t pos = 0;
    while ((pos = text.find(needle, pos)) != std::string::npos) {
        text.replace(pos, needle.size(), replacement);
        pos += replacement.size();
    }
    return text;
}

struct LegacyExpect {
    bool embedded = false;
    std::string filePath;
    std::vector<std::string> embeddedLines;
};

struct DirectiveSpec {
    int expectExit = 0;
    bool stdoutIgnore = true;
    std::string stdoutFile;
    std::vector<std::string> stderrSubstrings;
    std::string outputName;
    std::string gfaExpectPath;
    std::string gfaPreserveInput = "skip";
    std::string expectGfaHeader;
};

struct TestCase {
    std::string command;
    bool legacyMode = true;
    LegacyExpect legacy;
    DirectiveSpec directives;
};

struct TagValue {
    char type = '\0';
    std::string content;
};

struct GfaSegment {
    std::string name;
    std::string sequence;
    std::map<std::string, TagValue> tags;
};

struct GfaConnection {
    char recordType = '\0';
    std::string from;
    char fromOrient = '\0';
    std::string to;
    char toOrient = '\0';
    std::string payload;
    std::map<std::string, TagValue> tags;
};

struct GfaDocument {
    std::string headerVersion;
    std::multiset<std::string> nonTelomereRawLines;
    std::multiset<std::string> nonTelomereSegmentKeys;
    std::vector<GfaSegment> telomereSegments;
    std::vector<GfaConnection> telomereConnections;
};

struct GfaExpectRow {
    std::string path;
    std::string segment;
    std::string terminalRole;
    char pathOrient = '.';
    std::string nodeName;
    char segEdgeOrient = '\0';
    char connectorType = 'L';
    std::string connectorValue = "0M";
    int tlBp = 0;
};

bool parseTestFile(const std::string &path, TestCase &testCase, std::string &error) {
    std::ifstream stream(path);
    if (!stream) {
        error = "couldn't open test file";
        return false;
    }

    if (!std::getline(stream, testCase.command)) {
        error = "missing command line";
        return false;
    }
    testCase.command = stripCarriageReturn(testCase.command);

    std::vector<std::string> remainder;
    for (std::string line; std::getline(stream, line); )
        remainder.push_back(stripCarriageReturn(line));

    size_t firstNonEmpty = remainder.size();
    for (size_t i = 0; i < remainder.size(); ++i) {
        if (!trim(remainder[i]).empty()) {
            firstNonEmpty = i;
            break;
        }
    }

    if (firstNonEmpty == remainder.size()) {
        error = "missing expected output or directives";
        return false;
    }

    const std::string firstDirective = trim(remainder[firstNonEmpty]);
    const bool directiveMode = startsWith(firstDirective, "expect_") || startsWith(firstDirective, "gfa_");
    testCase.legacyMode = !directiveMode;

    if (!directiveMode) {
        const std::string legacyRef = trim(remainder[firstNonEmpty]);
        if (legacyRef == "embedded") {
            testCase.legacy.embedded = true;
            testCase.legacy.embeddedLines.assign(remainder.begin() + static_cast<long>(firstNonEmpty + 1), remainder.end());
        } else {
            testCase.legacy.embedded = false;
            testCase.legacy.filePath = legacyRef;
        }
        return true;
    }

    for (const std::string &raw : remainder) {
        const std::string line = trim(raw);
        if (line.empty() || startsWith(line, "#"))
            continue;

        std::istringstream iss(line);
        std::string key;
        iss >> key;
        std::string value;
        std::getline(iss, value);
        value = trim(value);

        if (key == "expect_exit") {
            testCase.directives.expectExit = std::stoi(value);
        } else if (key == "expect_stdout") {
            if (value == "ignore") {
                testCase.directives.stdoutIgnore = true;
            } else if (!value.empty()) {
                testCase.directives.stdoutIgnore = false;
                testCase.directives.stdoutFile = value;
            } else {
                error = "expect_stdout requires a value";
                return false;
            }
        } else if (key == "expect_stderr_substr") {
            testCase.directives.stderrSubstrings.push_back(value);
        } else if (key == "expect_output_name") {
            testCase.directives.outputName = value;
        } else if (key == "expect_gfa_header") {
            testCase.directives.expectGfaHeader = value;
        } else if (key == "gfa_expect") {
            testCase.directives.gfaExpectPath = value;
        } else if (key == "gfa_preserve_input") {
            testCase.directives.gfaPreserveInput = value;
        } else {
            error = "unknown directive: " + key;
            return false;
        }
    }

    return true;
}

bool compareLegacyExpected(std::istream &expectedStream, const fs::path &actualPath, const std::string &inputFile) {
    std::ifstream actOutput(actualPath);
    if (!actOutput) {
        printFAIL(inputFile.c_str(), "couldn't open captured stdout");
        return false;
    }

    std::string line;
    std::getline(expectedStream, line);
    if (line == "+++Summary+++: ") {
        std::getline(actOutput, line);
        std::set<std::string> expSummary, actSummary;
        while (!actOutput.eof()) {
            std::getline(actOutput, line);
            actSummary.insert(line);
        }
        while (!expectedStream.eof()) {
            std::getline(expectedStream, line);
            expSummary.insert(line);
        }

        std::set<std::string> additions, missings;
        for (const auto &entry : expSummary) {
            if (actSummary.count(entry) == 0)
                missings.insert(entry);
        }
        for (const auto &entry : actSummary) {
            if (expSummary.count(entry) == 0)
                additions.insert(entry);
        }

        if (!additions.empty() || !missings.empty()) {
            printFAIL(inputFile.c_str(), "expected output did not match actual output");
            std::cout << "additions:" << std::endl;
            for (const auto &addition : additions)
                std::cout << addition << std::endl;
            std::cout << "missing:" << std::endl;
            for (const auto &missing : missings)
                std::cout << missing << std::endl;
            return false;
        }
        return true;
    }

    std::vector<std::pair<std::string, std::string>> diffs;
    std::string actualLine;
    std::getline(actOutput, actualLine);
    if (actualLine != line)
        diffs.emplace_back(actualLine, line);

    std::string expectedLine;
    while (!actOutput.eof() || !expectedStream.eof()) {
        std::getline(actOutput, actualLine);
        std::getline(expectedStream, expectedLine);
        if (actualLine != expectedLine)
            diffs.emplace_back(actualLine, expectedLine);
    }

    if (!diffs.empty()) {
        printFAIL(inputFile.c_str(), "expected output did not match actual output");
        for (const auto &pair : diffs) {
            std::cout << "    expected: " << pair.second.c_str() << std::endl
                      << "      actual: " << pair.first.c_str() << std::endl;
        }
        return false;
    }

    return true;
}

std::map<std::string, TagValue> parseTags(const std::vector<std::string> &fields, size_t startIndex) {
    std::map<std::string, TagValue> tags;
    for (size_t i = startIndex; i < fields.size(); ++i) {
        std::vector<std::string> parts = split(fields[i], ':');
        if (parts.size() < 3 || parts[0].size() < 2)
            continue;
        TagValue value;
        value.type = parts[1].empty() ? '\0' : parts[1][0];
        value.content = joinTail(parts, 2, ':');
        tags[parts[0].substr(0, 2)] = value;
    }
    return tags;
}

std::string makeSegmentKey(const GfaSegment &segment) {
    std::ostringstream key;
    key << segment.name << '\t' << segment.sequence;
    for (const auto &tag : segment.tags)
        key << '\t' << tag.first << ':' << tag.second.type << ':' << tag.second.content;
    return key.str();
}

bool isTelomereName(const std::string &name) {
    return startsWith(name, "telomere_");
}

GfaDocument parseGfa(const fs::path &path) {
    GfaDocument doc;
    std::ifstream stream(path);
    if (!stream)
        throw std::runtime_error("couldn't open GFA file: " + path.string());

    std::string line;
    while (std::getline(stream, line)) {
        line = stripCarriageReturn(line);
        if (line.empty())
            continue;

        const std::vector<std::string> fields = split(line, '\t');
        if (fields.empty())
            continue;

        const std::string &recordType = fields[0];
        if (recordType == "H") {
            for (size_t i = 1; i < fields.size(); ++i) {
                if (startsWith(fields[i], "VN:Z:"))
                    doc.headerVersion = fields[i].substr(5);
            }
            doc.nonTelomereRawLines.insert(line);
            continue;
        }

        if (recordType == "S" && fields.size() >= 3) {
            const bool gfa2Segment = (!doc.headerVersion.empty() && doc.headerVersion[0] == '2');
            const std::string &name = fields[1];
            GfaSegment segment;
            segment.name = name;
            if (gfa2Segment && fields.size() >= 4) {
                segment.sequence = fields[3];
                segment.tags = parseTags(fields, 4);
            } else {
                segment.sequence = fields[2];
                segment.tags = parseTags(fields, 3);
            }

            if (isTelomereName(name)) {
                doc.telomereSegments.push_back(segment);
            } else {
                doc.nonTelomereRawLines.insert(line);
                doc.nonTelomereSegmentKeys.insert(makeSegmentKey(segment));
            }
            continue;
        }

        if ((recordType == "L" || recordType == "J") && fields.size() >= 6) {
            const bool telomereConnection = isTelomereName(fields[1]) || isTelomereName(fields[3]);
            GfaConnection connection;
            connection.recordType = recordType[0];
            connection.from = fields[1];
            connection.fromOrient = fields[2].empty() ? '\0' : fields[2][0];
            connection.to = fields[3];
            connection.toOrient = fields[4].empty() ? '\0' : fields[4][0];
            connection.payload = fields[5];
            connection.tags = parseTags(fields, 6);
            if (telomereConnection)
                doc.telomereConnections.push_back(connection);
            else
                doc.nonTelomereRawLines.insert(line);
            continue;
        }

        doc.nonTelomereRawLines.insert(line);
    }

    return doc;
}

std::vector<GfaExpectRow> parseGfaExpectations(const fs::path &path) {
    std::ifstream stream(path);
    if (!stream)
        throw std::runtime_error("couldn't open GFA expectation file: " + path.string());

    std::string line;
    std::unordered_map<std::string, size_t> header;
    std::vector<GfaExpectRow> rows;

    while (std::getline(stream, line)) {
        line = stripCarriageReturn(line);
        const std::string trimmed = trim(line);
        if (trimmed.empty() || startsWith(trimmed, "#"))
            continue;

        const std::vector<std::string> fields = split(line, '\t');
        if (header.empty()) {
            for (size_t i = 0; i < fields.size(); ++i)
                header[fields[i]] = i;
            continue;
        }

        auto getField = [&](const std::string &name) -> std::string {
            const auto it = header.find(name);
            if (it == header.end() || it->second >= fields.size())
                return "";
            return fields[it->second];
        };

        GfaExpectRow row;
        row.path = getField("path");
        row.segment = getField("segment");
        row.terminalRole = getField("terminal_role");
        const std::string orient = getField("path_orient");
        row.pathOrient = orient.empty() ? '.' : orient[0];
        row.nodeName = getField("node_name");
        const std::string segOrient = getField("seg_edge_orient");
        row.segEdgeOrient = segOrient.empty() ? '\0' : segOrient[0];
        const std::string connectorType = getField("connector_type");
        if (!connectorType.empty())
            row.connectorType = connectorType[0];
        const std::string connectorValue = getField("connector_value");
        if (!connectorValue.empty())
            row.connectorValue = connectorValue;
        const std::string tl = getField("tl_bp");
        row.tlBp = tl.empty() ? 0 : std::stoi(tl);
        rows.push_back(row);
    }

    return rows;
}

bool compareMultisetEquality(const std::multiset<std::string> &lhs,
                             const std::multiset<std::string> &rhs,
                             const std::string &inputFile,
                             const char *message) {
    if (lhs == rhs)
        return true;

    printFAIL(inputFile.c_str(), message);
    std::map<std::string, int> lhsCounts, rhsCounts;
    for (const auto &entry : lhs)
        lhsCounts[entry]++;
    for (const auto &entry : rhs)
        rhsCounts[entry]++;

    std::cout << "input-only:" << std::endl;
    for (const auto &[entry, count] : lhsCounts) {
        const int diff = count - rhsCounts[entry];
        for (int i = 0; i < diff; ++i)
            std::cout << "    " << entry << std::endl;
    }
    std::cout << "output-only:" << std::endl;
    for (const auto &[entry, count] : rhsCounts) {
        const int diff = count - lhsCounts[entry];
        for (int i = 0; i < diff; ++i)
            std::cout << "    " << entry << std::endl;
    }
    return false;
}

bool compareMultisetSubset(const std::multiset<std::string> &subset,
                           const std::multiset<std::string> &superset,
                           const std::string &inputFile,
                           const char *message) {
    std::map<std::string, int> subsetCounts, supersetCounts;
    for (const auto &entry : subset)
        subsetCounts[entry]++;
    for (const auto &entry : superset)
        supersetCounts[entry]++;

    std::vector<std::string> missing;
    for (const auto &[entry, count] : subsetCounts) {
        for (int i = supersetCounts[entry]; i < count; ++i)
            missing.push_back(entry);
    }

    if (missing.empty())
        return true;

    printFAIL(inputFile.c_str(), message);
    for (const auto &entry : missing)
        std::cout << "    missing in output: " << entry << std::endl;
    return false;
}

bool checkGfaPreservation(const DirectiveSpec &spec,
                          const GfaDocument &inputDoc,
                          const GfaDocument &outputDoc,
                          const std::string &inputFile) {
    if (spec.gfaPreserveInput == "skip" || spec.gfaPreserveInput.empty())
        return true;
    if (spec.gfaPreserveInput == "strict") {
        return compareMultisetEquality(inputDoc.nonTelomereRawLines,
                                       outputDoc.nonTelomereRawLines,
                                       inputFile,
                                       "non-telomere GFA lines changed");
    }
    if (spec.gfaPreserveInput == "subset") {
        return compareMultisetSubset(inputDoc.nonTelomereRawLines,
                                     outputDoc.nonTelomereRawLines,
                                     inputFile,
                                     "input GFA lines missing from annotated output");
    }
    if (spec.gfaPreserveInput == "segments_only") {
        return compareMultisetEquality(inputDoc.nonTelomereSegmentKeys,
                                       outputDoc.nonTelomereSegmentKeys,
                                       inputFile,
                                       "non-telomere GFA segments changed");
    }

    printFAIL(inputFile.c_str(), "unknown gfa_preserve_input mode", spec.gfaPreserveInput.c_str());
    return false;
}

std::string expectedNodeName(const GfaExpectRow &row) {
    if (row.pathOrient == '.')
        return "telomere_" + row.segment + "_" + row.terminalRole;
    return "telomere_" + row.segment + std::string(1, row.pathOrient) + "_" + row.terminalRole;
}

bool checkGfaExpectations(const GfaDocument &outputDoc,
                          const std::vector<GfaExpectRow> &expectedRows,
                          const std::string &inputFile) {
    std::multiset<std::string> actualNames;
    for (const auto &segment : outputDoc.telomereSegments)
        actualNames.insert(segment.name);

    std::multiset<std::string> expectedNames;
    for (const auto &row : expectedRows)
        expectedNames.insert(row.nodeName);

    if (actualNames != expectedNames) {
        printFAIL(inputFile.c_str(), "telomere node set mismatch");
        std::cout << "expected nodes:" << std::endl;
        for (const auto &name : expectedNames)
            std::cout << "    " << name << std::endl;
        std::cout << "actual nodes:" << std::endl;
        for (const auto &name : actualNames)
            std::cout << "    " << name << std::endl;
        return false;
    }

    if (outputDoc.telomereConnections.size() != expectedRows.size()) {
        printFAIL(inputFile.c_str(), "unexpected number of telomere connectors");
        std::cout << "    expected: " << expectedRows.size() << std::endl;
        std::cout << "    actual: " << outputDoc.telomereConnections.size() << std::endl;
        return false;
    }

    for (const auto &row : expectedRows) {
        const std::string expectedName = expectedNodeName(row);
        if (row.nodeName != expectedName) {
            printFAIL(inputFile.c_str(), "invalid expected node name", row.nodeName.c_str(), expectedName.c_str());
            return false;
        }

        std::vector<GfaSegment> matchingSegments;
        for (const auto &segment : outputDoc.telomereSegments) {
            if (segment.name == row.nodeName)
                matchingSegments.push_back(segment);
        }
        if (matchingSegments.size() != 1) {
            printFAIL(inputFile.c_str(), "expected exactly one telomere segment", row.nodeName.c_str());
            return false;
        }

        const auto &segment = matchingSegments.front();
        auto findTag = [&](const std::string &name) -> const TagValue* {
            const auto it = segment.tags.find(name);
            return it == segment.tags.end() ? nullptr : &it->second;
        };

        const TagValue *ln = findTag("LN");
        const TagValue *rc = findTag("RC");
        const TagValue *tl = findTag("TL");
        if (!ln || ln->content != "6" || !rc || rc->content != "6000" ||
            !tl || tl->content != std::to_string(row.tlBp)) {
            printFAIL(inputFile.c_str(), "telomere segment tags mismatch", row.nodeName.c_str());
            return false;
        }

        std::vector<GfaConnection> matchingConnections;
        for (const auto &connection : outputDoc.telomereConnections) {
            if (connection.from == row.nodeName)
                matchingConnections.push_back(connection);
        }
        if (matchingConnections.size() != 1) {
            printFAIL(inputFile.c_str(), "expected exactly one telomere connector", row.nodeName.c_str());
            return false;
        }

        const auto &connection = matchingConnections.front();
        if (connection.recordType != row.connectorType ||
            connection.fromOrient != '+' || connection.to != row.segment ||
            connection.toOrient != row.segEdgeOrient || connection.payload != row.connectorValue) {
            printFAIL(inputFile.c_str(), "telomere connector mismatch", row.nodeName.c_str());
            return false;
        }
        const auto rcTag = connection.tags.find("RC");
        if (rcTag == connection.tags.end() || rcTag->second.content != "0") {
            printFAIL(inputFile.c_str(), "telomere connector RC tag mismatch", row.nodeName.c_str());
            return false;
        }
    }

    return true;
}

std::string extractInputPath(const std::string &command) {
    const std::vector<std::string> tokens = split(command, ' ');
    for (size_t i = 0; i < tokens.size(); ++i) {
        if (tokens[i] == "-f" || tokens[i] == "--input-sequence") {
            if (i + 1 < tokens.size())
                return tokens[i + 1];
        }
    }
    for (const std::string &token : tokens) {
        if (!token.empty() && token[0] != '-')
            return token;
    }
    return "";
}

bool compareDirectiveStdout(const DirectiveSpec &spec,
                            const fs::path &capturedStdout,
                            const std::string &inputFile) {
    if (spec.stdoutIgnore)
        return true;

    std::ifstream expected(spec.stdoutFile);
    if (!expected) {
        printFAIL(inputFile.c_str(), "couldn't open expected stdout file", spec.stdoutFile.c_str());
        return false;
    }
    return compareLegacyExpected(expected, capturedStdout, inputFile);
}

bool runDirectiveAssertions(const TestCase &testCase,
                            const std::string &inputFile,
                            const std::string &command,
                            const fs::path &stdoutPath,
                            const fs::path &stderrPath,
                            const fs::path &outDir,
                            int actualExit) {
    const DirectiveSpec &spec = testCase.directives;
    bool ok = true;

    if (actualExit != spec.expectExit) {
        printFAIL(inputFile.c_str(), "unexpected exit code");
        std::cout << "    expected: " << spec.expectExit << std::endl;
        std::cout << "    actual: " << actualExit << std::endl;
        const std::string stderrText = readFile(stderrPath);
        if (!stderrText.empty())
            std::cout << stderrText << std::endl;
        return false;
    }

    ok = compareDirectiveStdout(spec, stdoutPath, inputFile) && ok;

    const std::string stderrText = readFile(stderrPath);
    for (const std::string &needle : spec.stderrSubstrings) {
        if (stderrText.find(needle) == std::string::npos) {
            printFAIL(inputFile.c_str(), "stderr substring missing", needle.c_str());
            ok = false;
        }
    }

    const std::string inputGfaPath = extractInputPath(command);
    const bool needsGfa = !spec.gfaExpectPath.empty() || !spec.expectGfaHeader.empty() ||
                          !spec.outputName.empty() || spec.gfaPreserveInput != "skip";
    if (!needsGfa)
        return ok;

    const std::string derivedOutput = fs::path(inputGfaPath).filename().string() + ".telo.annotated.gfa";
    const fs::path outputPath = outDir / (spec.outputName.empty() ? derivedOutput : spec.outputName);
    if (!fs::exists(outputPath)) {
        printFAIL(inputFile.c_str(), "expected GFA output missing", outputPath.string().c_str());
        return false;
    }

    try {
        const GfaDocument outputDoc = parseGfa(outputPath);
        if (!spec.expectGfaHeader.empty() && outputDoc.headerVersion != spec.expectGfaHeader) {
            printFAIL(inputFile.c_str(), "unexpected GFA header version", outputDoc.headerVersion.c_str());
            ok = false;
        }

        if (!inputGfaPath.empty() && spec.gfaPreserveInput != "skip") {
            const GfaDocument inputDoc = parseGfa(inputGfaPath);
            ok = checkGfaPreservation(spec, inputDoc, outputDoc, inputFile) && ok;
        }

        if (!spec.gfaExpectPath.empty()) {
            const std::vector<GfaExpectRow> expectedRows = parseGfaExpectations(spec.gfaExpectPath);
            ok = checkGfaExpectations(outputDoc, expectedRows, inputFile) && ok;
        }
    } catch (const std::exception &e) {
        printFAIL(inputFile.c_str(), "GFA assertion error", e.what());
        return false;
    }

    return ok;
}

} // namespace

int main(int argc, char **argv) {
    if (argc == 1) {
        std::cout << "teloscope-validate <path to test folder and/or files>" << std::endl;
        exit(EXIT_SUCCESS);
    }

    int opt;
    while ((opt = getopt(argc, argv, "c")) != -1) {
        switch (opt) {
        case 'c':
            printCommand = true;
            break;
        }
    }

    std::set<std::string> inputFiles;
    for (int i = 1; i < argc; ++i)
        get_recursive(argv[i], inputFiles);

    const std::string exePath = getExePath(argv[0]);

    for (const auto &inputFile : inputFiles) {
#ifdef _WIN32
        if (inputFile.find(".gz") != std::string::npos)
            continue;
#endif

        TestCase testCase;
        std::string parseError;
        if (!parseTestFile(inputFile, testCase, parseError)) {
            printFAIL(inputFile.c_str(), parseError.c_str());
            continue;
        }

        const auto uniqueSeed = std::chrono::steady_clock::now().time_since_epoch().count();
        const fs::path baseDir = fs::temp_directory_path() /
                                 ("teloscope_validate_" + sanitizeName(inputFile) + "_" + std::to_string(uniqueSeed));
        const fs::path outDir = baseDir / "out";
        const fs::path stdoutPath = baseDir / "stdout.txt";
        const fs::path stderrPath = baseDir / "stderr.txt";
        fs::create_directories(baseDir);

        const bool needsOutDir =
            testCase.command.find("%OUTDIR%") != std::string::npos ||
            !testCase.directives.outputName.empty();
        std::string command = replaceAll(testCase.command, "%OUTDIR%", outDir.string());
        if (needsOutDir)
            fs::create_directories(outDir);

#ifdef _WIN32
        const std::string systemCommand =
            "\"\"" + exePath + "\"" + " " + command + " > " + stdoutPath.string() + " 2>" + stderrPath.string() + "\"";
#else
        const std::string systemCommand =
            shellQuote(exePath) + " " + command + " > " + shellQuote(stdoutPath.string()) +
            " 2>" + shellQuote(stderrPath.string());
#endif

        if (printCommand)
            std::cout << systemCommand << std::endl;

        const int status = system(systemCommand.c_str());
        const int actualExit = normalizeExitCode(status);

        bool ok = true;
        if (testCase.legacyMode) {
            if (actualExit != EXIT_SUCCESS) {
                printFAIL(inputFile.c_str(), "runtime error");
                const std::string stderrText = readFile(stderrPath);
                if (!stderrText.empty())
                    std::cout << stderrText << std::endl;
                ok = false;
            } else if (testCase.legacy.embedded) {
                std::stringstream embeddedExpected;
                for (const std::string &line : testCase.legacy.embeddedLines)
                    embeddedExpected << line << '\n';
                ok = compareLegacyExpected(embeddedExpected, stdoutPath, inputFile);
            } else {
                std::ifstream expected(testCase.legacy.filePath);
                if (!expected) {
                    printFAIL(inputFile.c_str(), "couldn't open expected output", testCase.legacy.filePath.c_str());
                    ok = false;
                } else {
                    ok = compareLegacyExpected(expected, stdoutPath, inputFile);
                }
            }
        } else {
            ok = runDirectiveAssertions(testCase, inputFile, command, stdoutPath, stderrPath, outDir, actualExit);
        }

        fs::remove_all(baseDir);

        if (ok)
            printPASS(inputFile.c_str());
    }

    exit(pass ? EXIT_SUCCESS : EXIT_FAILURE);
}
