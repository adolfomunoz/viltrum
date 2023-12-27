#pragma once
#include "regions-integrator-sequential.h"
#include "../integrator-per-bin-parallel.h"

namespace viltrum {

/**
 * @brief This integrator is not intended to be directly used by the user. Is a particular piece
 * common to all integrators that are based on Newton-Cotes regions. 
 */
class RegionsIntegratorParallelBins {
public:
	template<typename Bins, std::size_t DIMBINS, typename SeqRegions, typename Float, std::size_t DIM, typename Logger>
	void integrate_regions(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const SeqRegions& seq_regions, const Range<Float,DIM>& range, Logger& logger) const {
        
        auto dummy_function = 
            [] (const std::array<Float,DIM>& x) { return typename std::decay_t<decltype(seq_regions.front())>::value_type(); };

        integrator_per_bin_parallel(regions_integrator_sequential().integrator_wrapper(seq_regions)).
            integrate(bins,bin_resolution,dummy_function,range,logger);
    }
};


/**
 * @brief This integrator is not intended to be directly used by the user. Is a particular piece
 * common to all integrators that are based on Newton-Cotes regions. The integrand is therefore ignored
 * 
 * @tparam SeqRegions Any STL-like linear sequence of regions
 */
auto regions_integrator_parallel_bins() {
    return RegionsIntegratorParallelBins();
}

}