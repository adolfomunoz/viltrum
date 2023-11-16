#include "../../src/monte-carlo/random-sequence.h"
#include "../../src/combination/concat.h"
#include <iostream>
#include <list>



using namespace viltrum;

int main() {
    auto l = concat(std::list{5.0,8.0},random_sequence(RangeInfinite(std::vector{-2.0},std::vector{2.0})));
    int i = 0;
    for (double d : l) {
        std::cout<<d<<" ";
        if (++i>=10) break;
    } 
    std::cout<<std::endl;
}  