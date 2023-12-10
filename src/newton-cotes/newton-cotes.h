#pragma once
#include "region.h"
#include "integrator-regions.h"
#include "rules.h"

namespace viltrum {

template<typename Rule>
class NewtonCotes {
    Rule rule;
    bool parallel;

public:
    NewtonCotes(const Rule& r, bool parallel = false) : rule(r),parallel(parallel) {} 

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {
        auto logger_region = logger_step(logger,"integrand evaluation");
        logger_region.log_progress(0,1);
        auto r = region(f,rule, range.min(), range.max());
        logger_region.log_progress(1,1);

        auto logger_integration = logger_step(logger, "subrange integration");
        integrate_regions(bins,bin_resolution,std::vector<std::decay_t<decltype(r)>>{r},range,parallel,logger_integration);
	}
};

template<typename R>
NewtonCotes<R> newton_cotes(const R& rule, bool parallel = false){
    return NewtonCotes<R>(rule, parallel);
} 



}