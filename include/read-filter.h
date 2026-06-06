#ifndef READ_FILTER_H
#define READ_FILTER_H

#include <memory>
#include <string>

struct UserInputTeloscope;
class Teloscope;

class ReadTelomereFilter {
    std::unique_ptr<Teloscope> teloscope;

public:
    explicit ReadTelomereFilter(const UserInputTeloscope &input);
    ~ReadTelomereFilter();
    bool matches(std::string sequence);
};

#endif /* READ_FILTER_H */
