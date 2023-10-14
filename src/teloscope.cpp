#include <istream>
#include <fstream>
#include <sstream>

#include <parallel_hashmap/phmap.h>

#include "log.h"
#include "global.h"
#include "uid-generator.h"

#include "bed.h"
#include "struct.h"
#include "functions.h"

#include "gfa-lines.h"
#include "gfa.h"
#include "sak.h"

#include "stream-obj.h"

#include "input-agp.h"
#include "input-filters.h"
#include "input-gfa.h"

#include "teloscope.h"

BedCoordinates bedCoords;

std::vector<uint64_t> maxSum(std::vector<bool> seq, uint32_t windowSize, int step){
    std::vector<uint64_t> aver;
    uint64_t pSize = seq.size();
    aver.reserve(pSize - 5);

    if (pSize < windowSize){
        return aver;
    }

    uint32_t max_sum = 0; // check this
    for (uint32_t i = 0; i < windowSize; ++i){
        if (seq[i])
            max_sum += 1;
    }
    uint32_t window_sum = max_sum;
    aver.push_back(window_sum);

    for (uint64_t i = windowSize; i < pSize; i += step){
        window_sum += seq[i] - seq[i - windowSize]; // rolling average -> do this 
        max_sum = std::max(max_sum, window_sum);
        aver.push_back(window_sum);
    }
    
    return aver; // float? change for count
}

void findTelomeres(std::string pHeader, std::string &P, UserInput userInput){

    std::string S = userInput.telomerePattern;
    uint16_t sSize = S.size();
    uint64_t pSize = P.size();
    int windowSize = userInput.windowSize; // change int for.. less variables
    int step = userInput.step; 
        std::vector<bool> telo_location (pSize, false);
    std::vector<uint64_t> aver;
    std::vector<uint32_t> pos;
        

	for (uint64_t i = 0;i < pSize - sSize + 1; ++i){
		if (S == P.substr(i, sSize)) {
            telo_location[i] = true;
            pos.push_back(i);
            i = i + sSize - 1;
		}
    }

    aver = maxSum(telo_location, windowSize, step);

    // for (uint64_t value : pos) {
    //     std::cout << pHeader << "  ";
    //     std::cout << value << "  " << value + windowSize - 1 << "  ";
    //     std::cout << aver[value] << std::endl;
    // }

    // for (uint64_t value : pos) {
    //     std::cout << pHeader << "  ";
    //     std::cout << value <<"  "<< value + sSize - 1 << "  ";
    //     std::cout << std::endl;
    // }

    for (uint64_t value : pos) {
        bedCoords.pushCoordinates(pHeader, value, value + sSize - 1);
    }

    std::string bedFileName = "./output/" + userInput.nameInput + "_telomere_locations.bed";

    std::ofstream bedFile(bedFileName);
    for (unsigned int i = 0; i < bedCoords.size(); ++i) {
        bedFile << bedCoords.getSeqHeader(i) << "\t" << bedCoords.getBegin(i) << "\t" << bedCoords.getEnd(i) << std::endl; // replace to /n 
    }
    bedFile.close();

}