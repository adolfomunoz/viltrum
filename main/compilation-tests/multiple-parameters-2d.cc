#include "../../viltrum.h"
#include <iostream>
#include <iomanip>

using namespace viltrum;

int main() {
    const std::size_t bins = 10;
    const unsigned long samples = 1000000;

    auto f =[] (float x, float y,float z) -> double {
        if (((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))<0.25) return 1.0;
        else return 0.0;
    };

    //With floats, due to numerica stability, for a very large number of samples, monte_carlo fails to converge?
    {
        LoggerProgress logger("Simple");
        std::cout<<integrate(monte_carlo(samples*bins*bins),f,range_primary<3>(),logger)<<std::endl;  
    }

    {
        LoggerProgress logger("Monte-Carlo");
        std::vector<std::vector<float>> sol(bins,std::vector<float>(bins,0.0f));
        integrate(monte_carlo(bins*bins*samples),sol,f,range_primary<3>(),logger);
        for (const auto& vv : sol) {
            for (float v : vv)
                std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }
} 