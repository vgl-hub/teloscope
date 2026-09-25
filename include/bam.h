#ifndef BAM_H
#define BAM_H

#include <iosfwd>

struct UserInputTeloscope;
struct ReadTlStats;

void readBamReads(const UserInputTeloscope &userInput, std::ostream &subset,
                  std::ostream &bed, ReadTlStats &stats);

#endif /* BAM_H */
