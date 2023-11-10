#include "../../viltrum.h"
#include <iostream>
#include <iomanip>

using namespace viltrum;

int main() {
    const std::size_t bins = 24;
    const unsigned long samples = 1000000;

    auto f =[] (const std::array<float,2>& x) -> double {
       if ((x[0]+x[1])<1) return 1.0;        
       else return 0.0;
    };

    //With floats, due to numerica stability, for a very large number of samples, monte_carlo fails to converge?
    {
        LoggerProgress logger("Simple");
        std::cout<<integrate(monte_carlo(samples*bins),f,range_primary<2>(),logger)<<std::endl;  
    }

    {
        LoggerProgress logger("Monte-Carlo");
        std::vector<float> sol(bins,0.0f);
        integrate(monte_carlo(bins*samples),sol,f,range_primary<2>(),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }
    {
        LoggerProgress logger("Per bin");
        std::vector<float> sol(bins,0.0f); 
        integrate(integrator_per_bin(monte_carlo(samples)),sol,f,range_primary<2>(),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }
    {
        LoggerProgress logger("Parallel"); //I'm suprised this works with a RNG shared among threads
        std::vector<float> sol(bins,0.0f); 
        integrate(integrator_per_bin_parallel(monte_carlo(samples)),sol,f,range_primary<2>(),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }
    {
        LoggerProgress logger("Trapezoids");
        std::vector<float> sol(bins,0.0f); 
        integrate(newton_cotes(steps<bins/2>(trapezoidal)),sol,f,range_primary<2>(),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }
    {
        LoggerProgress logger("Perbin trapezoids");
        std::vector<float> sol(bins,0.0f); 
        integrate(integrator_per_bin(newton_cotes(steps<16>(trapezoidal))),sol,f,range_primary<2>(),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }

} 