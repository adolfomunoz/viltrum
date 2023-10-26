#include "../viltrum.h"
#include <iostream>

using namespace viltrum;

int main() {
    const std::size_t bins = 16;
    const unsigned long samples = 409600;
    auto f =[] (const std::array<float,2>& x) {
        if ((x[0]+x[1])<1) return 1.0f;
        else return 0.0f;
    };
    {
        LoggerProgress logger("Simple");
        std::cout<<integrate(monte_carlo(samples),f,range_primary<2>(),logger)<<std::endl;  
    }

    {
        LoggerProgress logger("Monte-Carlo");
        std::vector<float> sol(bins,0.0f);
        integrate(monte_carlo(bins*samples),sol,f,range_primary<2>(),logger);
        for (float v : sol) std::cout<<v<<" ";
        std::cout<<std::endl;
    }
    {
        LoggerProgress logger("Per bin");
        std::vector<float> sol(bins,0.0f); 
        integrate(integrator_per_bin(monte_carlo(samples)),sol,f,range_primary<2>(),logger);
        for (float v : sol) std::cout<<v<<" ";
        std::cout<<std::endl;
    }

} 