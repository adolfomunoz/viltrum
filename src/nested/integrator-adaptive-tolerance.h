#pragma once
#include "../newton-cotes/region.h"
#include "../newton-cotes/regions-integrator-sequential.h"
#include "nested.h"
#include "error-heuristic.h"

namespace viltrum {

template<typename Rule, typename ErrorHeuristic, typename = std::enable_if_t<is_nested<Rule>::value>>
class IntegratorAdaptiveTolerance {
    Rule rule;
    ErrorHeuristic error_heuristic;
    float tolerance;

    template<typename R, typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate_region(const R& r, Float& integrated_volume, Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {
            auto [error,dimension] = error_heuristic(r);
            if (error < tolerance) {
                LoggerNull logger_region;
                regions_integrator_sequential().integrate_regions(bins,bin_resolution,std::array<R,1>{r},range,logger_region);
                integrated_volume += r.range().volume();
                logger.log_progress(integrated_volume,range.volume());
            } else {
                auto subregions = r.split(f,dimension);
                for (const auto& subr : subregions) this->integrate_region(subr,integrated_volume,bins,bin_resolution,f,range,logger);
            }
        }
public:
    IntegratorAdaptiveTolerance(const Rule& r, const ErrorHeuristic& er, float tol) : 
        rule(r),error_heuristic(er), tolerance(tol) {} 

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {
            Float integrated_volume(0);
            auto r = region(f,rule, range.min(), range.max());
            integrate_region(r,integrated_volume,bins,bin_resolution,f,range,logger);
        }
};

template<typename Rule, typename ErrorHeuristic>
auto integrator_adaptive_tolerance(const Rule& rule, const ErrorHeuristic& error_heuristic, float tolerance = 1.e-3) {
    return IntegratorAdaptiveTolerance<Rule,ErrorHeuristic>(rule,error_heuristic, tolerance);
}

template<typename Rule>
auto integrator_adaptive_tolerance(const Rule& rule, float tolerance = 1.e-3) {
    return integrator_adaptive_tolerance(rule,error_heuristic_default(error_metric_absolute()),tolerance);
}

template<typename Rule>
auto integrator_adaptive_tolerance(const Rule& rule, double tolerance = 1.e-3) {
    return integrator_adaptive_tolerance(rule,error_heuristic_default(error_metric_absolute()),tolerance);
}


}