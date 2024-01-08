#pragma once
#include "../newton-cotes/region.h"
#include "../newton-cotes/regions-integrator-sequential.h"
#include "../newton-cotes/regions-integrator-parallel-regions.h"
#include "regions-generator-adaptive-heap.h"
#include "../newton-cotes/integrator-region-based.h"
#include "nested.h"
#include "error-heuristic.h"

namespace viltrum {

template<typename R, typename EH>
auto integrator_adaptive_iterations(const R& rule, const EH& error_heuristic, std::size_t iterations) {
    return integrator_region_based(regions_generator_adaptive_heap(rule,error_heuristic,iterations),regions_integrator_sequential());
}

template<typename R, typename EH>
auto integrator_adaptive_iterations_parallel(const R& rule, const EH& error_heuristic, std::size_t iterations, std::size_t nmutexes=16) {
    return integrator_region_based(regions_generator_adaptive_heap(rule,error_heuristic,iterations),regions_integrator_parallel_regions(nmutexes));
}

template<typename Rule>
auto integrator_adaptive_iterations(const Rule& rule, std::size_t iterations) {
    return integrator_adaptive_iterations(rule,error_heuristic_default(error_metric_absolute()),iterations);
}

template<typename Rule>
auto integrator_adaptive_iterations_parallel(const Rule& rule, std::size_t iterations) {
    return integrator_adaptive_iterations_parallel(rule,error_heuristic_default(error_metric_absolute()),iterations);
}



}