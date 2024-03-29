#include "../../viltrum.h"
#include <iostream>
#include <iomanip>

using namespace viltrum;

int main() {
    const std::size_t bins = 24;
    const unsigned long samples = 10000;

    auto f =[] (const auto& seq) -> double {
        auto it = seq.begin();
        double n = 0.0;
        while (true) {
            float x = *it; ++it;
            float y = *it; ++it;  
            if ((x*x+y*y)<1) return n;
            n+=1.0;
        } 
    };

/*
    //With floats, due to numerica stability, for a very large number of samples, monte_carlo fails to converge?
    {
        LoggerProgress logger("Simple");
        std::cout<<integrate(monte_carlo(samples*bins),f,range_infinite(-1.0,-1.0,1.0,1.0),logger)<<std::endl;  
    }

    {
        LoggerProgress logger("Monte-Carlo");
        std::vector<float> sol(bins,0.0f);
        integrate(monte_carlo(bins*samples),sol,f,range_infinite(-1.0,-1.0,1.0,1.0),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }
    
    {
        LoggerProgress logger("Per bin");
        std::vector<float> sol(bins,0.0f); 
        integrate(integrator_per_bin(monte_carlo(samples)),sol,f,range_infinite(-1.0,-1.0,1.0,1.0),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }
*/

    {
        LoggerProgress logger("Parallel"); //I'm suprised this works with a RNG shared among threads
        std::vector<float> sol(bins,0.0f); 
        integrate(integrator_per_bin_parallel(monte_carlo(100*samples)),sol,f,range_infinite(-1.0,-1.0,1.0,1.0),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }

/*
    {
        LoggerProgress logger("Fubini"); 
        std::vector<float> sol(bins,0.0f); 
        integrate(integrator_fubini<2>(integrator_newton_cotes(steps<bins/2>(trapezoidal)),monte_carlo(samples)),sol,f,range_infinite(-1.0,-1.0,1.0,1.0),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }

    {
        LoggerProgress logger("Fubini adaptive"); 
        std::vector<float> sol(bins,0.0f); 
        integrate(integrator_fubini<4>(integrator_adaptive_tolerance(nested(simpson,trapezoidal),1.e-4f),monte_carlo(16)),sol,f,range_infinite(-1.0,-1.0,1.0,1.0),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }
*/

    {
        LoggerProgress logger("PPA IS");
        std::vector<float> sol(bins,0.0f); 
        integrate(integrator_fubini({0,1,4,5},integrator_adaptive_variance_reduction_parallel(nested(simpson,trapezoidal),128,rr_integral_region(),cv_fixed_weight(0.0),samples),monte_carlo(10)),sol,f,range_infinite(-1.0,-1.0,1.0,1.0),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }

    {
        LoggerProgress logger("PPA CV");
        std::vector<float> sol(bins,0.0f); 
        integrate(integrator_fubini({0,1,4,5},integrator_adaptive_variance_reduction_parallel(nested(simpson,trapezoidal),128,rr_uniform_region(),cv_fixed_weight(1.0),samples),monte_carlo(10)),sol,f,range_infinite(-1.0,-1.0,1.0,1.0),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }

    {
        LoggerProgress logger("New Fubini");
        std::vector<float> sol(bins,0.0f); 
        integrate(integrator_adaptive_fubini_variance_reduction_parallel<4>(
                    nested(simpson,trapezoidal),error_heuristic_size(error_metric_absolute(),1.e-5),128,10,
                    rr_error_region(),cv_optimize_weight(),region_sampling_uniform(),samples*10),
                sol,f,range_infinite(-1.0,-1.0,1.0,1.0),logger);
        for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
        std::cout<<std::endl;
    }
    
} 