#pragma once
#include "error-metric.h"

namespace viltrum {

/**
 * A error heuristic uses a region and obtains a 2 elements tuple with the error (arithmetic) and the dimension (std::size_t) in which it occurs.
 **/

template<typename ErrorMetric>
class error_heuristic_default {
    ErrorMetric error_metric;
public:
    error_heuristic_default(const ErrorMetric& em) : error_metric(em) {}
    template<typename R>
    auto operator()(const R& region) const {
        return region.max_error_dimension(error_metric);
    }
};


}