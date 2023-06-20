#include <iostream>
#include <string>
#include <iostream>
#include <vector>




int main(){
	std::string S = "TTAGGG";
    uint8_t sSize = S.size();
    int pSize = 54;
    int step = 1;
    int windowSize = 15;
    int thresh = 1;
    std::vector<bool> telo_location (pSize, false);
    std::vector<uint64_t> aver;
    std::vector<int> pos;
    std::vector<bool> tt (pSize, false);

    pos.push_back(0);

    for (uint64_t i = 0; i < pSize; ++i){
        if (aver[i] >= thresh)
        tt[i] = true;
    }

    int mid = pSize/2;
    for (uint64_t i = 0; i < pSize; ++i){
        if (tt[i] == true){
            if (i > mid){
            int endP = tt[i + windowSize - 1];
            pos.push_back(endP);
            }
            if (i < mid){
            int startP = tt[i + sSize - 1];
            pos.push_back(startP);
            }
        }
    }

    pos.push_back(pSize);

    for (int values : pos){
        std::cout << values << ", ";
    }

}