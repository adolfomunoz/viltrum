#pragma once
#include "integrator-adaptive-variance-reduction.h"
#include "integrator-adaptive-fubini-variance-reduction.h"

namespace viltrum {

auto integrator_crespo2021(std::size_t iterations, std::size_t spp, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    
    return integrator_region_based(
        regions_generator_adaptive_heap(
            nested(simpson,trapezoidal),
            error_heuristic_size(error_metric_relative(),1.e-5),
            iterations),
        regions_integrator_parallel_variance_reduction(
            rr_uniform_region(), 
            cv_optimize_weight(), 
            region_sampling_uniform(), 
            std::mt19937(seed),
            spp,
            nmutexes)
    ); 
}

template<std::size_t DIM>
auto integrator_crespo2021_infinite(std::size_t iterations, std::size_t mc_samples, std::size_t spp, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    
    return integrator_region_based(
        regions_generator_fubini<DIM>(
            regions_generator_adaptive_heap(
                nested(simpson,trapezoidal),
                error_heuristic_size(error_metric_relative(),1.e-5),
                iterations
            ),
            monte_carlo(mc_samples, 2*seed + 1)
        ),
        regions_integrator_parallel_variance_reduction(
            rr_uniform_region(), 
            cv_optimize_weight(), 
            region_sampling_uniform(), 
            std::mt19937(seed),
            spp,
            nmutexes)
    ); 
}


}