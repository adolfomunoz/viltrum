#pragma once
#include "../newton-cotes/integrator-region-based.h"
#include "../combination/regions-generator-fubini.h"
#include "../nested/regions-generator-adaptive-heap.h"
#include "regions-integrator-debug-count-regions.h"
#include "../monte-carlo/monte-carlo.h"



namespace viltrum {

template<std::size_t N, typename R, typename EH, typename RNG>
auto integrator_adaptive_debug_count_regions(const R& rule, const EH& error_heuristic, std::size_t iterations, unsigned long mc_samples, RNG&& rng) {
    RNG rn = rng;
    return integrator_region_based(
        regions_generator_fubini<N>(regions_generator_adaptive_heap(rule,error_heuristic,iterations),monte_carlo(rn,mc_samples)),
        regions_integrator_debug_count_regions());
}

template<std::size_t N, typename R, typename EH>
auto integrator_adaptive_debug_count_regions(const R& rule, const EH& error_heuristic, std::size_t iterations, unsigned long mc_samples, std::size_t seed = std::random_device()()) {
    return integrator_adaptive_debug_count_regions<N>(rule,error_heuristic,iterations,mc_samples,std::mt19937(seed));
}

template<typename R, typename EH, typename RNG>
auto integrator_adaptive_debug_count_regions(const R& rule, const EH& error_heuristic, std::size_t iterations) {
    return integrator_region_based(
        regions_generator_adaptive_heap(rule,error_heuristic,iterations),
        regions_integrator_debug_count_regions());
}


}