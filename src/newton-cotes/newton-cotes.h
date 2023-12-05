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
        
        auto r = region(f,rule, range.min(), range.max());
        integrate_regions(bins,bin_resolution,std::vector{r},range,parallel,logger);
	}
};

template<typename R>
NewtonCotes<R> newton_cotes(const R& rule, bool parallel = false){
    return NewtonCotes<R>(rule, parallel);
} 

}