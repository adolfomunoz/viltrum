#pragma once
#include "range.h"
#include "multidimensional-range.h"
#include "integrate.h"

namespace viltrum {

template<typename Integrator>
class IntegratorPerBin {
    Integrator bin_integrator;

public:
	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range) const {
        double factor(1);   
        for (std::size_t i = 0; i < DIMBINS; ++i) factor*=bin_resolution[i];
        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(bin_resolution[i]);
        for (auto pos : multidimensional_range(bin_resolution)) {
            Range<Float,DIM> subrange = range;
            for (std::size_t i=0;i<DIMBINS;++i)
                subrange = subrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);
            bins(pos) = factor*viltrum::integrate(bin_integrator,f,subrange);
        }
	}

    IntegratorPerBin(Integrator&& pi) : 
	    bin_integrator(std::forward<Integrator>(pi)) { }
    IntegratorPerBin(const Integrator& pi) : 
	    bin_integrator(pi) { }

};

template<typename Integrator>
auto integrator_per_bin(Integrator&& i) {
    return IntegratorPerBin<std::decay_t<Integrator>>(std::forward<Integrator>(i));
}

}