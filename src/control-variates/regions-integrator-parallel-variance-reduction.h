#pragma once
#include <random>
#include <execution>
#include "../range.h"
#include "mutexed-tensor-vector.h"
#include "../foreach.h"
#include "region-russian-roulette.h"
#include "weight-strategy.h"



namespace viltrum {

template<typename RR, typename CV, typename RNG>
class RegionsIntegratorParallelVarianceReduction {
    RR rr;
    CV cv;
    mutable RNG rng;
    unsigned long samples;
    std::size_t nmutexes;

public:
    RegionsIntegratorParallelVarianceReduction(RR&& rr, CV&& cv, RNG&& r, unsigned long s,std::size_t n = 16) 
        : rr(rr), cv(cv), rng(r), samples(s), nmutexes(n) {}
	template<typename Bins, std::size_t DIMBINS, typename SeqRegions, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate_regions(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const SeqRegions& seq_regions, const F& f, const Range<Float,DIM>& range, Logger& logger) const {
        
        using Reg = typename SeqRegions::value_type;

        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(bin_resolution[i]);
        std::size_t factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) factor*=bin_resolution[i];
        MutexedTensorVector<const Reg*,DIMBINS> perbin(bin_resolution,nmutexes);
        std::size_t progress = 0; std::size_t final_progress = seq_regions.size();
        auto logger_bins = logger_step(logger, "region bin stratification");
        logger_bins.log_progress(progress,final_progress);
        std::for_each(std::execution::par_unseq, seq_regions.begin(),seq_regions.end(),[&] (const auto& r) {
            for (auto pos : pixels_in_region(r,bin_resolution,range)) {
                perbin.push_at(pos,&r);
            } 
        });  
        logger_bins.log_progress(final_progress,final_progress);
        auto logger_control_variates = logger_step(logger,"residual and variance reduction");
        for_each(parallel,multidimensional_range(bin_resolution),
            [&] (const std::array<std::size_t, DIMBINS>& pos) {
                using Sample = decltype(f(std::declval<std::array<Float,DIM>>()));
                Range<Float,DIM> binrange = range;
                for (std::size_t i=0;i<DIMBINS;++i)
                    binrange = binrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);
                std::vector<std::tuple<const Reg*,Range<Float,DIM>>> regions_ranges(perbin[pos].size(),std::tuple((const Reg*)nullptr,range_primary<DIM,Float>()));

                //This is the control variate and the region_bin_ranges
                Sample approximation(0);
                for (std::size_t i = 0; i<regions_ranges.size(); ++i) {
                    regions_ranges[i] = std::tuple(perbin[pos][i],binrange.intersection(perbin[pos][i]->range()));
                    if (!std::get<1>(regions_ranges[i]).empty())   
                        approximation += double(factor)*std::get<0>(regions_ranges[i])->integral_subrange(std::get<1>(regions_ranges[i]));
	            }
                //These are the samples for the residual, accumulated into accumulator
                auto accumulator = cv.accumulator(Sample(0));
                auto roulette = rr.russian_roulette(regions_ranges);
                for (unsigned long i=0;i<samples;++i) {
                    auto [chosen,rrfactor] = roulette.choose(rng);
                    const auto& [r, region_bin_range]  = regions_ranges[chosen];
                    std::array<Float,DIM> sample;
	                for (std::size_t i=0;i<DIM;++i) {
		                std::uniform_real_distribution<Float> dis(region_bin_range.min(i),region_bin_range.max(i));
		                sample[i] = dis(rng);
                    }
                    accumulator.push(
                        f(sample)*double(factor)*rrfactor*region_bin_range.volume(),
                        r->approximation_at(sample)*double(factor)*rrfactor*region_bin_range.volume()
                    );
                }
                bins(pos) = accumulator.integral(approximation);    
        }   ,logger_control_variates);
    }
};

template<typename RR, typename CV, typename RNG>
auto regions_integrator_parallel_variance_reduction(RR&& rr, CV&& cv, RNG&& rng, unsigned long samples, std::size_t nmutexes = 16, std::enable_if_t<!std::is_integral_v<RNG>,int> dummy = 0) {
    return RegionsIntegratorParallelVarianceReduction<RR,CV,RNG>(std::forward<RR>(rr), std::forward<CV>(cv), std::forward<RNG>(rng),samples,nmutexes);
}

template<typename RR, typename CV>
auto regions_integrator_parallel_variance_reduction(RR&& rr, CV&& cv, unsigned long samples, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return regions_integrator_parallel_variance_reduction(std::forward<RR>(rr), std::forward<CV>(cv), std::mt19937(seed),samples,nmutexes);
}

template<typename RR>
auto regions_integrator_parallel_variance_reduction(RR&& rr, double alpha, unsigned long samples, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return regions_integrator_parallel_variance_reduction(std::forward<RR>(rr), cv_fixed_weight(alpha), std::mt19937(seed),samples,nmutexes);
}

}