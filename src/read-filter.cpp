#include <algorithm>
#include <limits>

#include "main.h"
#include "functions.h"
#include "teloscope.h"
#include "read-filter.h"

namespace {

UserInputTeloscope makeReadFilterInput(const UserInputTeloscope &input) {
    UserInputTeloscope readInput = input;
    readInput.minBlockLen = 42; // the subset floor is fixed at 7 repeats; -l only shapes the measure scan

    // max/2 keeps the full read terminal without overflowing scanSegment's doubled limit.
    readInput.terminalLimit = std::numeric_limits<uint32_t>::max() / 2;
    readInput.terminalTolerance = std::numeric_limits<uint32_t>::max() / 2;
    readInput.ultraFastMode = true;
    readInput.outFasta = false;
    readInput.outWinRepeats = false;
    readInput.outGC = false;
    readInput.outEntropy = false;
    readInput.outMatches = false;
    readInput.outPlotReport = false;
    readInput.manualCuration = false;
    return readInput;
}

} // namespace

// measure scan: a read is scanned like a contig end, -l keeps its assembly default, only tip allowance and tile size change
UserInputTeloscope makeReadTlInput(const UserInputTeloscope &input) {
    UserInputTeloscope readInput = input;
    if (!readInput.terminalToleranceSet) {
        readInput.terminalTolerance = 300; // read tip allowance
    }
    if (!readInput.terminalLimitSet) {
        readInput.terminalLimit = 2000; // read tile size
    }
    // the first tile must cover the tip allowance, or the anchor search zone gets clamped short
    readInput.terminalLimit = std::max(readInput.terminalLimit, readInput.terminalTolerance);
    readInput.ultraFastMode = true;
    readInput.outFasta = false;
    readInput.outWinRepeats = false;
    readInput.outGC = false;
    readInput.outEntropy = false;
    readInput.outMatches = false;
    readInput.outPlotReport = false;
    readInput.manualCuration = false;
    return readInput;
}

ReadTelomereFilter::ReadTelomereFilter(const UserInputTeloscope &input)
    : teloscope(std::make_unique<Teloscope>(makeReadFilterInput(input))) {}

ReadTelomereFilter::~ReadTelomereFilter() = default;

bool ReadTelomereFilter::matches(std::string &sequence) {
    if (!sequence.empty() && sequence.back() == '\r') {
        sequence.pop_back();
    }
    unmaskSequence(sequence);

    SegmentData segmentData = teloscope->scanSegment(sequence, 0, true, true, true, true);
    return !segmentData.terminalBlocks.empty();
}

ReadTelomereScanner::ReadTelomereScanner(const UserInputTeloscope &input)
    : teloscope(std::make_unique<Teloscope>(makeReadTlInput(input))) {}

ReadTelomereScanner::~ReadTelomereScanner() = default;

std::vector<TelomereBlock> ReadTelomereScanner::scan(std::string &sequence) {
    if (!sequence.empty() && sequence.back() == '\r') {
        sequence.pop_back();
    }
    unmaskSequence(sequence);

    SegmentData segmentData = teloscope->scanSegment(sequence, 0, true, true, true, true);
    return std::move(segmentData.terminalBlocks);
}
