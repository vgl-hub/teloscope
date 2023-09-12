// Tasks: jack
// replace P for genome or similar 
// replace S for pattern or similar TELOMERE_SEQUENCE
// replace aver for average
// replace pos for positions
// S has to be dinamically defined (from common k-mers) or from array
// UserInput userInput not used!!!!
// pHeader for p.header
// L67: int windowSize = 500; shouldn't be hardcoded but an argument of teloscope

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

std::vector<uint64_t> maxSum(std::vector<bool> seq, uint32_t windowSize, int step){
    std::vector<uint64_t> aver;
    uint64_t pSize = seq.size();
    aver.reserve(pSize - 5);

    if (pSize < windowSize){
        return aver;
    }

    uint32_t max_sum = 0;
    for (uint32_t i = 0; i < windowSize; ++i){
        if (seq[i])
            max_sum += 1;
    }
    uint32_t window_sum = max_sum;
    aver.push_back(window_sum);

    for (uint64_t i = windowSize; i < pSize; i += step){
        window_sum += seq[i] - seq[i - windowSize];
        max_sum = std::max(max_sum, window_sum);
        aver.push_back(window_sum);
    }
    
    return aver;
}

void findTelomeres(std::string pHeader, std::string &P, UserInput userInput){

	std::string S = "TTAGGG";
    uint16_t sSize = S.size();
    uint64_t pSize = P.size();
    int step = 1;
    int windowSize = 500;
    ///int thresh = 1;
    std::vector<bool> telo_location (pSize, false);
    std::vector<uint64_t> aver;
    std::vector<uint32_t> pos;
    ///std::vector<bool> tt (pSize, false);

    

	for (uint64_t i = 0;i < pSize - sSize + 1; ++i){
		if (S == P.substr(i,6)) {
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

    for (uint64_t value : pos) {
        std::cout << pHeader << "  ";
        std::cout << value <<"  "<< value + sSize - 1 << "  ";
        std::cout << std::endl;
    }

    // jack: BED file
    std::ofstream bedFile("telomere_locations.bed"); // no hardcode!!!! use {genome}.fa > {genome}.bed

    for (uint64_t value : pos) {
        bedFile << pHeader << "\t" << value << "\t" << value + sSize << std::endl;
    }

    bedFile.close();
}



