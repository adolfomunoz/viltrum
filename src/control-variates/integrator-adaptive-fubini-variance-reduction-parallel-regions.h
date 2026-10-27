#pragma once
#include "regions-integrator-parallel-regions-variance-reduction.h"
#include "../combination/regions-generator-fubini.h"
#include "../monte-carlo/monte-carlo.h"


namespace viltrum {



template<std::size_t N, typename Stratification, typename RS, typename R, typename EH, typename RNG>
auto integrator_adaptive_fubini_variance_reduction_parallel_regions(const R& rule, const EH& error_heuristic, std::size_t iterations, unsigned long mc_samples, Stratification&& stratification, RS&& region_sampler, RNG&& rng, unsigned long spp, double alpha = 1.0, std::size_t nmutexes = 1024) {
    RNG rn = rng;
    return integrator_region_based(
        regions_generator_fubini<N>(regions_generator_adaptive_parallel(rule,error_heuristic,iterations),monte_carlo(rn,mc_samples)),
        regions_integrator_parallel_regions_variance_reduction(std::forward<Stratification>(stratification),std::forward<RS>(region_sampler), rn,spp,alpha,nmutexes));
}

template<std::size_t N, typename Stratification, typename RS, typename R, typename EH>
auto integrator_adaptive_fubini_variance_reduction_parallel_regions(const R& rule, const EH& error_heuristic, std::size_t iterations, unsigned long mc_samples, Stratification&& stratification, RS&& region_sampler, unsigned long spp, double alpha = 1.0, std::size_t seed = std::random_device()(), std::size_t nmutexes = 1024) {
    return integrator_region_based(
        regions_generator_fubini<N>(regions_generator_adaptive_parallel(rule,error_heuristic,iterations),monte_carlo(std::mt19937(seed),mc_samples)),
        regions_integrator_parallel_regions_variance_reduction(std::forward<Stratification>(stratification),std::forward<RS>(region_sampler), std::mt19937(seed+1),spp,alpha,nmutexes));
}



}