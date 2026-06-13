#include <limits>

#include "main.h"
#include "functions.h"
#include "teloscope.h"
#include "read-filter.h"

namespace {

UserInputTeloscope makeReadFilterInput(const UserInputTeloscope &input) {
    UserInputTeloscope readInput = input;
    if (!readInput.minBlockLenSet) {
        readInput.minBlockLen = 42; // ~7 telomeric repeats
    }

    // max/2 keeps the full read terminal without overflowing scanSegment's doubled limit.
    readInput.terminalLimit = std::numeric_limits<uint32_t>::max() / 2;
    readInput.ultraFastMode = true;
    readInput.outFasta = false;
    readInput.outWinRepeats = false;
    readInput.outGC = false;
    readInput.outEntropy = false;
    readInput.outMatches = false;
    readInput.outITS = false;
    readInput.outPlotReport = false;
    readInput.manualCuration = false;
    return readInput;
}

} // namespace

ReadTelomereFilter::ReadTelomereFilter(const UserInputTeloscope &input)
    : teloscope(std::make_unique<Teloscope>(makeReadFilterInput(input))) {}

ReadTelomereFilter::~ReadTelomereFilter() = default;

bool ReadTelomereFilter::matches(std::string sequence) {
    if (!sequence.empty() && sequence.back() == '\r') {
        sequence.pop_back();
    }
    unmaskSequence(sequence);

    SegmentData segmentData = teloscope->scanSegment(sequence, 0, true);
    return !segmentData.terminalBlocks.empty();
}
