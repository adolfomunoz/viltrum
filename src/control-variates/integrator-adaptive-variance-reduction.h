#pragma once
#include "../newton-cotes/region.h"
#include "../nested/regions-generator-adaptive-heap.h"
#include "../newton-cotes/integrator-region-based.h"
#include "../nested/nested.h"
#include "../nested/error-heuristic.h"
#include "regions-integrator-parallel-variance-reduction.h"

namespace viltrum {

template<typename RR, typename R, typename EH, typename RNG>
auto integrator_adaptive_variance_reduction_parallel(const R& rule, const EH& error_heuristic, std::size_t iterations, RNG&& rng, unsigned long spp, std::size_t nmutexes = 16) {
    return integrator_region_based(regions_generator_adaptive_heap(rule,error_heuristic,iterations),regions_integrator_parallel_variance_reduction<RR>(std::forward<RNG>(rng),spp,nmutexes));
}

template<typename RR, typename R, typename EH>
auto integrator_adaptive_variance_reduction_parallel(const R& rule, const EH& error_heuristic, std::size_t iterations, unsigned long spp, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return integrator_region_based(regions_generator_adaptive_heap(rule,error_heuristic,iterations),regions_integrator_parallel_variance_reduction<RR>(spp,seed,nmutexes));
}

template<typename RR, typename R>
auto integrator_adaptive_variance_reduction_parallel(const R& rule, std::size_t iterations, unsigned long spp, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return integrator_adaptive_variance_reduction_parallel<RR>(rule,error_heuristic_default(error_metric_absolute()),iterations,spp,seed,nmutexes);
}


}