#include "../../src/monte-carlo/random-sequence.h"
#include <iostream>



using namespace viltrum;

int main() {
    auto l = random_sequence(RangeInfinite(std::vector{-2.0},std::vector{2.0}));
    int i = 0;
    for (double d : l) {
        std::cout<<d<<" ";
        if (++i>=10) break;
    } 
    std::cout<<std::endl;
}  