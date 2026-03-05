#pragma once
#include <random>
#include <execution>
#include "../range.h"
#include "../foreach.h"
#include "../combination/fubini.h"
#include "../monte-carlo/monte-carlo.h"
#include "../monte-carlo/monte-carlo-per-bin.h"
#include <array>

namespace viltrum {

class RegionsIntegratorDebugCountRegions {

public:
    RegionsIntegratorDebugCountRegions() {}

	template<typename Bins, std::size_t DIMBINS, typename SeqRegions, typename F, typename IntegrationRange, typename Logger>
	void integrate_regions(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const SeqRegions& seq_regions, const F& f, const IntegrationRange& range, Logger& logger) const {
        
        using Reg = typename SeqRegions::value_type;
        using RegRange = std::decay_t<decltype(std::declval<Reg>().range())>;
        constexpr std::size_t DIM = RegRange::size;
        using Float = typename RegRange::value_type;
        auto ranges = range_split_at<DIM>(range);
        auto range_first = std::get<0>(ranges); auto range_rest = std::get<1>(ranges);
        // vv This would be in C++20
        // auto [range_first,range_rest] = range_split_at<DIM>(range); 

        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(bin_resolution[i]);
        std::size_t factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) factor*=bin_resolution[i];
        std::size_t progress = 0; std::size_t final_progress = seq_regions.size();
        auto logger_bins = logger_step(logger, "region counting");
        logger_bins.log_progress(progress,final_progress);
        auto f_regdim = function_split_and_integrate_at<DIM>(f,monte_carlo(1),range_rest);
        using Sample = std::decay_t<decltype(f_regdim(std::declval<std::array<Float,DIM>>()))>;
        

        std::for_each(std::execution::seq, seq_regions.begin(),seq_regions.end(),[&] (const auto& r) {
            for (auto pos : pixels_in_region(r,bin_resolution,range_first)) {
                bins(pos) += Sample(1);
            } 
        }); 
        //Maybe we need to debug that "pixels_in_region" is working correctly, checking if the
        // intersection between the bin and the region is never empty.
    }
};

auto regions_integrator_debug_count_regions() {
    return RegionsIntegratorDebugCountRegions();
}

}