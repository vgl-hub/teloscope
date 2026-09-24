#ifndef READ_FILTER_H
#define READ_FILTER_H

#include <memory>
#include <string>
#include <vector>

struct UserInputTeloscope;
class Teloscope;
struct TelomereBlock;

// the measure scan's effective parameters (tip allowance, tile size, tile/tolerance floor), also written as provenance
UserInputTeloscope makeReadTlInput(const UserInputTeloscope &input);

class ReadTelomereFilter {
    std::unique_ptr<Teloscope> teloscope;

public:
    explicit ReadTelomereFilter(const UserInputTeloscope &input);
    ~ReadTelomereFilter();
    bool matches(std::string sequence);
};

// read TL mode: same engine as ReadTelomereFilter, but returns the terminal blocks themselves
class ReadTelomereScanner {
    std::unique_ptr<Teloscope> teloscope;

public:
    explicit ReadTelomereScanner(const UserInputTeloscope &input);
    ~ReadTelomereScanner();
    std::vector<TelomereBlock> scan(std::string sequence);
};

#endif /* READ_FILTER_H */
