#include "../../src/monte-carlo/random-sequence.h"
#include <iostream>



using namespace viltrum;

int main() {
    auto l = random_sequence(0.0,1.0);
    int i = 0;
    for (double d : l) {
        std::cout<<d<<" ";
        if (++i>=10) break;
    } 
    std::cout<<std::endl;
}  