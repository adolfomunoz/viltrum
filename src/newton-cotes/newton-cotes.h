#pragma once
#include "region.h"
#include "rules.h"

namespace viltrum {

template<typename Rule>
class NewtonCotes {
    Rule rule;

public:
    NewtonCotes(const Rule& r) : rule(r) {} 

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {
        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(bin_resolution[i]);
        std::size_t factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) factor*=bin_resolution[i];
        auto r = region(f,rule, range.min(), range.max());
        std::size_t i = 0;

        for (auto pos : multidimensional_range(bin_resolution)) {
            logger.log_progress(i++,factor);
            auto subrange = range;
            for (std::size_t i=0;i<DIMBINS;++i)
                subrange = subrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);
            bins(pos) += double(factor)*r.integral_subrange(subrange);
        }
        logger.log_progress(factor,factor);
	}
};

template<typename R>
NewtonCotes<R> newton_cotes(const R& rule){
    return NewtonCotes<R>(rule);
} 

}