#pragma once
#include <random>
#include "../range.h"

namespace viltrum {

template<typename RNG>
class RegionsIntegratorParallelControlVariatesLowMemory {
    mutable RNG rng;
    unsigned long samples;
    std::size_t nmutexes;

public:
    RegionsIntegratorParallelControlVariatesLowMemory(RNG&& r, unsigned long s,std::size_t n = 16) 
        : rng(r), samples(s), nmutexes(n) {}
	template<typename Bins, std::size_t DIMBINS, typename SeqRegions, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate_regions(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const SeqRegions& seq_regions, const F& f, const Range<Float,DIM>& range, Logger& logger) const {
        
        using Reg = typename SeqRegions::value_type;

        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(bin_resolution[i]);
        std::size_t factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) factor*=bin_resolution[i];
        std::size_t progress = 0; std::size_t final_progress = seq_regions.size();

        auto logger_control_variates = logger_step(logger,"control variates");
        for_each(parallel,multidimensional_range(bin_resolution),
            [&] (const std::array<std::size_t, DIMBINS>& pos) {
                Range<Float,DIM> binrange = range;
                for (std::size_t i=0;i<DIMBINS;++i)
                    binrange = binrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);
                
                std::vector<std::tuple<const Reg*,Range<Float,DIM>>> regions_ranges;
                std::for_each(seq_regions.begin(),seq_regions.end(),[&] (const auto& r) {
                    auto region_bin_range = binrange.intersection(r.range());
                    if (!region_bin_range.empty()) {
                        regions_ranges.push_back(std::tuple(&r,region_bin_range));
                    }
                });

                //This is for the control variate 
                for (const auto& [r,region_bin_range] : regions_ranges) {
                    bins(pos) += double(factor)*r->integral_subrange(region_bin_range);
                }
                //This is the MonteCarlo residual, RR among regions
                std::uniform_int_distribution<std::size_t> rr(0,regions_ranges.size()-1);
                for (unsigned long i=0;i<samples;++i) {
                    std::size_t chosen = rr(rng);
                    const auto& [r,region_bin_range] = regions_ranges[chosen];
                    std::array<Float,DIM> sample;
                    for (std::size_t i=0;i<DIM;++i) {
                        std::uniform_real_distribution<Float> dis(region_bin_range.min(i),region_bin_range.max(i));
                        sample[i] = dis(rng);
                    }
                    bins(pos) += (f(sample)-r->approximation_at(sample))*double(factor)*region_bin_range.volume()*double(regions_ranges.size())/double(samples);
                }
        }   ,logger_control_variates);
    }
};

template<typename RNG>
auto regions_integrator_parallel_control_variates_low_memory(RNG&& rng, unsigned long samples, std::size_t nmutexes = 16, std::enable_if_t<!std::is_integral_v<RNG>,int> dummy = 0) {
    return RegionsIntegratorParallelControlVariatesLowMemory<RNG>(std::forward<RNG>(rng),samples,nmutexes);
}

auto regions_integrator_parallel_control_variates_low_memory(unsigned long samples, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return regions_integrator_parallel_control_variates_low_memory(std::mt19937(seed),samples,nmutexes);
}

}