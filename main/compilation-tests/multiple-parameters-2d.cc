#include "../../viltrum.h"
#include <iostream>
#include <iomanip>

using namespace viltrum;

int main() {
    const std::size_t bins = 10;
    const unsigned long samples = 10000;

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

/*
    {
        LoggerProgress logger("Monte-Carlo per-bin parallel");
        std::vector<std::vector<float>> sol(bins,std::vector<float>(bins,0.0f));
        integrate(integrator_per_bin_parallel(monte_carlo(samples)),sol,f,range_primary<3>(),logger);
        for (const auto& vv : sol) {
            for (float v : vv)
                std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }

    {
        LoggerProgress logger("Control variates");
        std::vector<std::vector<float>> sol(bins,std::vector<float>(bins,0.0f));
        integrate(integrator_adaptive_control_variates_parallel_low_memory(nested(simpson,trapezoidal),128,samples),sol,f,range_primary<3>(),logger);
        for (const auto& vv : sol) {
            for (float v : vv)
                std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }

    {
        LoggerProgress logger("Adaptive variance reduction - full importance sampling");
        std::vector<std::vector<float>> sol(bins,std::vector<float>(bins,0.0f));
        integrate(integrator_adaptive_variance_reduction_parallel(nested(simpson,trapezoidal),128,rr_integral_region(),cv_fixed_weight(0.0),samples),sol,f,range_primary<3>(),logger);
        for (const auto& vv : sol) {
            for (float v : vv)
                std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }

    {
        LoggerProgress logger("Adaptive variance reduction - optimize cv-is weight");
        std::vector<std::vector<float>> sol(bins,std::vector<float>(bins,0.0f));
        integrate(integrator_adaptive_variance_reduction_parallel(nested(simpson,trapezoidal),128,rr_integral_region(),cv_optimize_weight(),samples),sol,f,range_primary<3>(),logger);
        for (const auto& vv : sol) {
            for (float v : vv)
                std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }

    {
        LoggerProgress logger("New Fubini");
        std::vector<std::vector<float>> sol(bins,std::vector<float>(bins,0.0f));
        integrate(integrator_adaptive_fubini_variance_reduction_parallel<2>(
                    nested(simpson,trapezoidal),error_heuristic_default(error_metric_absolute()),128,10,
                    rr_integral_region(),cv_optimize_weight(),region_sampling_uniform(),samples),
                sol,f,range_primary<3>(),logger);
        for (const auto& vv : sol) {
            for (float v : vv)
                std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }
*/
    {
        LoggerProgress logger("New Fubini Importance");
        std::vector<std::vector<float>> sol(bins,std::vector<float>(bins,0.0f));
        integrate(integrator_adaptive_fubini_variance_reduction<2>(
                    nested(simpson,trapezoidal),error_heuristic_default(error_metric_absolute()),128,10,
                    rr_integral_region(),cv_fixed_weight(0),region_sampling_importance(),samples),
                sol,f,range_primary<3>(),logger);
        for (const auto& vv : sol) {
            for (float v : vv)
                std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }
} 