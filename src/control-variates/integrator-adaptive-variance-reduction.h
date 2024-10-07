#pragma once
#include "../newton-cotes/region.h"
#include "../nested/regions-generator-adaptive-heap.h"
#include "../newton-cotes/integrator-region-based.h"
#include "../nested/nested.h"
#include "../nested/error-heuristic.h"
#include "regions-integrator-parallel-variance-reduction.h"
#include "regions-integrator-variance-reduction.h"

namespace viltrum {

template<typename RR, typename CV, typename RS, typename R, typename EH, typename RNG>
auto integrator_adaptive_variance_reduction_parallel(const R& rule, const EH& error_heuristic, std::size_t iterations, RR&& rr, CV&& cv, RS&& region_sampler, RNG&& rng, unsigned long spp, std::size_t nmutexes = 16) {
    return integrator_region_based(regions_generator_adaptive_heap(rule,error_heuristic,iterations),regions_integrator_parallel_variance_reduction(std::forward<RR>(rr),std::forward<CV>(cv), std::forward<RS>(region_sampler), std::forward<RNG>(rng),spp,nmutexes));
}

template<typename RR, typename CV, typename R, typename EH, typename RNG>
auto integrator_adaptive_variance_reduction_parallel(const R& rule, const EH& error_heuristic, std::size_t iterations, RR&& rr, CV&& cv, RNG&& rng, unsigned long spp, std::size_t nmutexes = 16) {
    return integrator_region_based(regions_generator_adaptive_heap(rule,error_heuristic,iterations),regions_integrator_parallel_variance_reduction(std::forward<RR>(rr),std::forward<CV>(cv), std::forward<RNG>(rng),spp,nmutexes));
}

template<typename RR, typename CV, typename R, typename EH>
auto integrator_adaptive_variance_reduction_parallel(const R& rule, const EH& error_heuristic, std::size_t iterations, RR&& rr, CV&& cv, unsigned long spp, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return integrator_region_based(regions_generator_adaptive_heap(rule,error_heuristic,iterations),regions_integrator_parallel_variance_reduction(std::forward<RR>(rr),std::forward<CV>(cv), spp,seed,nmutexes));
}

template<typename RR, typename CV, typename RS, typename R, typename EH>
auto integrator_adaptive_variance_reduction_parallel(const R& rule, const EH& error_heuristic, std::size_t iterations, RR&& rr, CV&& cv, RS&& region_sampler, unsigned long spp, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return integrator_region_based(regions_generator_adaptive_heap(rule,error_heuristic,iterations),regions_integrator_parallel_variance_reduction(std::forward<RR>(rr),std::forward<CV>(cv), std::forward<RS>(region_sampler), spp,seed,nmutexes));
}

template<typename RR, typename CV, typename RS, typename R>
auto integrator_adaptive_variance_reduction_parallel(const R& rule, std::size_t iterations, RR&& rr, CV&& cv, RS&& region_sampler, unsigned long spp, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return integrator_adaptive_variance_reduction_parallel(rule,error_heuristic_default(error_metric_absolute()),iterations,std::forward<RR>(rr),std::forward<CV>(cv), std::forward<RS>(region_sampler),spp,seed,nmutexes);
}



template<typename RR, typename CV, typename R>
auto integrator_adaptive_variance_reduction_parallel(const R& rule, std::size_t iterations, RR&& rr, CV&& cv, unsigned long spp, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return integrator_adaptive_variance_reduction_parallel(rule,error_heuristic_default(error_metric_absolute()),iterations,std::forward<RR>(rr),std::forward<CV>(cv), spp,seed,nmutexes);
}


}