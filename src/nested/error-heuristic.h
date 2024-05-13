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

template<typename ErrorMetric>
class error_heuristic_size {
    ErrorMetric error_metric;
    double size_weight; double min_size;
public:
    error_heuristic_size(const ErrorMetric& em, double sw = 1.e-5, double ms = 1.e-37) : 
            error_metric(em), size_weight(sw), min_size(ms) {}

    template<typename R>
    auto operator()(const R& region) const {
        auto max_err = region.error(0,error_metric);
        if ((region.range().max(0)-region.range().min(0))<min_size) max_err=0;
        max_err += size_weight*std::abs(region.range().max(0) - region.range().min(0)); 
        std::size_t max_dim = 0; 
        auto err = max_err;
		for (std::size_t d = 1; d<R::dimensions; ++d) {
			err = region.error(d,error_metric) + 
                    size_weight*std::abs(region.range().max(d) - region.range().min(d));

            if ((region.range().max(d)-region.range().min(d))<min_size) err=0;
			if (err>=max_err) {
				max_err = err; max_dim = d;
			}
		}
		return std::tuple(max_err,max_dim);
    } 
};



}