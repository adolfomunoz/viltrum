#include "../viltrum.h"
#include <iostream>




using namespace viltrum;

int main() {
    auto f =[] (const std::array<float,2>& x) {
        if ((x[0]+x[1])<1) return 1.0f;
        else return 0.0f;
    };
    std::cout<<integrate(monte_carlo(16),f,range_primary<2>())<<std::endl;  

} 