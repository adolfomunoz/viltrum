#pragma once
#include "region.h"
#include "regions-integrator-sequential.h"
#include "regions-integrator-parallel-bins.h"
#include "regions-generator-single.h"
#include "integrator-region-based.h"
#include "rules.h"

namespace viltrum {

template<typename R>
auto integrator_newton_cotes(const R& rule) {
    return integrator_region_based(regions_generator_single(rule),regions_integrator_sequential());
}

template<typename R>
auto integrator_newton_cotes_parallel(const R& rule) {
    return integrator_region_based(regions_generator_single(rule),regions_integrator_parallel_bins());
}



}