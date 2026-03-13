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
        if ((region.range().max(0)-region.range().min(0))<min_size || std::isnan(region.range().max(0)-region.range().min(0))) max_err=0;
        else max_err += size_weight*std::abs(region.range().max(0) - region.range().min(0)); 
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

template<typename ErrorMetricBins, typename ErrorMetricRest>
class error_heuristic_mixed {
    ErrorMetricBins error_metric_bins;
    ErrorMetricRest error_metric_rest;
    unsigned int dimension; //Dimension that sepparates the bins from the rest of the dimensions. The error metric for the bins is applied to this dimension, and the error metric for the rest is applied to the other dimensions.
    double bins_weight; //Relative weight of the error metric for the bins compared to the error metric for the rest of the dimensions. The error is computed as bins_weight*error_metric_bins + (1-bins_weight)*error_metric_rest.
    double size_threshold_bins; //Size threshold for the bins dimensions. Should be the inverse of the resolution (product for all dimensions). If the size of the region in this dimension is smaller than this threshold, the error extremelly increased.
    double size_threshold_rest; //Size threshold for the rest of the dimensions. Should be the inverse of the spp. If the size of the region in any of these dimensions is smaller than this threshold, the error extremelly increased.
    double error_increase_factor; //Factor by which the error is increased when the size of the region in any dimension is smaller than the corresponding threshold. The error is multiplied by this factor.
    double size_weight; //Relative weight of the size of the region compared to the error. The error is computed as error + size_weight*size per dimension.
public:
    error_heuristic_mixed(const ErrorMetricBins& em_bins, const ErrorMetricRest& em_rest, 
                    unsigned int dim = 2, double bins_w = 1.0, double size_weight = 1.e-3,
                    double size_bins = 1.0/1024.0, double size_rest = 1.0/16.0,
                    double error_increase_factor = 1.e4) : 
            error_metric_bins(em_bins), error_metric_rest(em_rest), dimension(dim), 
            bins_weight(bins_w), size_threshold_bins(size_bins), size_threshold_rest(size_rest), error_increase_factor(error_increase_factor), size_weight(size_weight) {}

    template<typename R>
    auto operator()(const R& region) const {
        double size_bins = 1.0; double size_rest = 1.0;
        for (std::size_t d = 0; d<(std::min(dimension,(unsigned int)(R::dimensions))); ++d) 
            size_bins *= std::abs(region.range().max(d)-region.range().min(d));
        for (std::size_t d = dimension; d<R::dimensions; ++d)
            size_rest *= std::abs(region.range().max(d)-region.range().min(d));

        double add_error_for_bins = error_increase_factor; 
        double add_error_for_rest = error_increase_factor;
        if (size_bins<size_threshold_bins) add_error_for_bins = 0.0;
        if (size_rest<size_threshold_rest) add_error_for_rest = 0.0;
        if (std::isnan(size_bins)) add_error_for_bins = 0.0;
        if (std::isnan(size_rest)) add_error_for_rest = 0.0;

        auto max_err = region.error(0,error_metric_bins)*bins_weight + 
            add_error_for_bins + size_weight*(region.range().max(0)-region.range().min(0));
        std::size_t max_dim = 0;
        auto err = max_err;
//        std::cerr<<std::setprecision(12)<<std::scientific<<"("<<size_bins<<","<<size_rest<<")  "<<"0 "<<err<<" | ";
		for (std::size_t d = 1; d<R::dimensions; ++d) {
            if (d<dimension) err = region.error(d,error_metric_bins)*bins_weight + add_error_for_bins  + size_weight*(region.range().max(d)-region.range().min(d));
            else err = region.error(d,error_metric_rest) + add_error_for_rest  + size_weight*(region.range().max(d)-region.range().min(d));
//            std::cerr<<d<<" "<<err<<" | ";
			if (err>=max_err) {
				max_err = err; max_dim = d;
			}
		}
//        std::cerr<<" -> "<<max_err<<" "<<max_dim<<std::endl;
		return std::tuple(max_err,max_dim);
    } 
};



}