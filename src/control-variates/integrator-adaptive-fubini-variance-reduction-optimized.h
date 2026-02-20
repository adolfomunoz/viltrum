#pragma once
#include "integrator-adaptive-variance-reduction-optimized.h"
#include "../combination/regions-generator-fubini.h"
#include "../monte-carlo/monte-carlo.h"

namespace viltrum {


template<std::size_t N, typename RR, typename CV, typename RS, typename R, typename EH, typename RNG>
auto integrator_adaptive_fubini_variance_reduction_parallel_optimized(const R& rule, const EH& error_heuristic, std::size_t iterations, unsigned long mc_samples, RR&& rr, CV&& cv, RS&& region_sampler, RNG&& rng, unsigned long spp, std::size_t nmutexes = 16) {
    RNG rn = rng;
    return integrator_region_based(
        regions_generator_fubini<N>(regions_generator_adaptive_heap(rule,error_heuristic,iterations),monte_carlo(rn,mc_samples)),
        regions_integrator_parallel_variance_reduction_optimized(std::forward<RR>(rr),std::forward<CV>(cv), std::forward<RS>(region_sampler), rn,spp,nmutexes));
}

template<std::size_t N, typename RR, typename CV, typename RS, typename R, typename EH>
auto integrator_adaptive_fubini_variance_reduction_parallel_optimized(const R& rule, const EH& error_heuristic, std::size_t iterations, unsigned long mc_samples, 
        RR&& rr, CV&& cv, RS&& region_sampler, unsigned long spp, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return integrator_region_based(
        regions_generator_fubini<N>(regions_generator_adaptive_heap(rule,error_heuristic,iterations),monte_carlo(std::mt19937(seed),mc_samples)),
        regions_integrator_parallel_variance_reduction_optimized(std::forward<RR>(rr),std::forward<CV>(cv), std::forward<RS>(region_sampler), std::mt19937(seed+1),spp,nmutexes));
}


}