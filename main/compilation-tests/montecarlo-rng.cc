#include <iostream>
#include <iomanip>
#include "../../viltrum.h"
#include <cmath>

class DysonSeries {
    float decay;
public:
    DysonSeries(float decay = 0.75f) : decay(decay) {}
    template<typename Seq>
    float operator()(const Seq& seq) const {
        float sum = 1.0f; 
        float range = 1.0f;
        auto it = seq.begin();
        while ((*it) < decay) { //Multiply by decay, divide by decay as probability, gets simplified.
            ++it;
            float x = range*(*it);
            float prob = 1.0/range;
            range = x; //The new x reduces the range
            sum += x/prob;
            ++it;
        }
        return sum;
    }
} integrand_infinite;

int main(int argc, char **argv) {
    auto range_infinite = viltrum::range_primary_infinite<float>();
    std::size_t nbins=10;
    unsigned long samples = 100000;

    for (int i = 0;i<(argc-1);++i) {
        if (std::string(argv[i])=="-bins") nbins = std::atoi(argv[++i]);
        else if (std::string(argv[i])=="-samples") samples = std::atoi(argv[++i]);
    } 

    {
        viltrum::LoggerProgress logger("Default");
        std::vector<float> sol_bins(nbins,0.0f);
        viltrum::integrate(viltrum::integrator_per_bin_parallel(viltrum::monte_carlo(samples,0)), sol_bins, integrand_infinite, range_infinite, logger);
        for (float x : sol_bins) std::cout<<x<<" ";
        std::cout<<std::endl;
    }

    {
        viltrum::LoggerProgress logger("Xoshiro128PlusPlus");
        std::vector<float> sol_bins(nbins,0.0f);
        viltrum::integrate(viltrum::integrator_per_bin_parallel(viltrum::monte_carlo(XoshiroCpp::Xoshiro128PlusPlus(0), samples)), sol_bins, integrand_infinite, range_infinite, logger);
        for (float x : sol_bins) std::cout<<x<<" ";
        std::cout<<std::endl;
    }  
    
    {
        viltrum::LoggerProgress logger("pcg32");
        std::vector<float> sol_bins(nbins,0.0f);
        viltrum::integrate(viltrum::integrator_per_bin_parallel(viltrum::monte_carlo(pcg32(0), samples)), sol_bins, integrand_infinite, range_infinite, logger);
        for (float x : sol_bins) std::cout<<x<<" ";
        std::cout<<std::endl;
    }  

	return 0;
}



