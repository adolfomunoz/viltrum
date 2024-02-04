#pragma once
#include "range.h"
#include "multidimensional-range.h"
#include "foreach.h"
#include "integrate.h"


namespace viltrum {

template<typename Integrator>
class IntegratorPerBinParallel {
    Integrator bin_integrator;

public:
	template<typename Bins, std::size_t DIMBINS, typename F, typename R, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const R& range, Logger& logger) const {
        std::size_t factor(1);   
        for (std::size_t i = 0; i < DIMBINS; ++i) factor*=bin_resolution[i];
        std::array<typename R::value_type,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/((typename R::value_type)(bin_resolution[i]));

        for_each(parallel,multidimensional_range(bin_resolution),
            [&] (const std::array<std::size_t,DIMBINS>& pos) {
                R subrange = range;
                for (std::size_t i=0;i<DIMBINS;++i)
                    subrange = subrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);
                bins(pos) = double(factor)*viltrum::integrate(bin_integrator,f,subrange);
            },logger);
	}

    IntegratorPerBinParallel(Integrator&& pi) : 
	    bin_integrator(std::forward<Integrator>(pi)) { }
    IntegratorPerBinParallel(const Integrator& pi) : 
	    bin_integrator(pi) { }

};

template<typename Integrator>
auto integrator_per_bin_parallel(Integrator&& i) {
    return IntegratorPerBinParallel<std::decay_t<Integrator>>(std::forward<Integrator>(i));
}

}