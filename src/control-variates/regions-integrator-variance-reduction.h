#pragma once
#include <random>
#include <execution>
#include "../range.h"
#include "mutexed-tensor-vector.h"
#include "../foreach.h"
#include "region-russian-roulette.h"
#include "weight-strategy.h"
#include "region-sampling.h"
#include "../combination/fubini.h"
#include "../monte-carlo/monte-carlo.h"

namespace viltrum {

template<typename RR, typename CV, typename RS, typename RNG>
class RegionsIntegratorVarianceReduction {
    RR rr;
    CV cv;
    RS region_sampler;
    mutable RNG rng;
    unsigned long samples;
    std::size_t nmutexes;

public:
    RegionsIntegratorVarianceReduction(RR&& rr, CV&& cv, RS&& rs, RNG&& r, unsigned long s,std::size_t n = 16) 
        : rr(rr), cv(cv), region_sampler(rs), rng(r), samples(s), nmutexes(n) {}
	template<typename Bins, std::size_t DIMBINS, typename SeqRegions, typename F, typename IntegrationRange, typename Logger>
	void integrate_regions(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const SeqRegions& seq_regions, const F& f, const IntegrationRange& range, Logger& logger) const {
        
        using Reg = typename SeqRegions::value_type;
        using RegRange = std::decay_t<decltype(std::declval<Reg>().range())>;
        constexpr std::size_t DIM = RegRange::size;
        using Float = typename RegRange::value_type;
        auto [range_first,range_rest] = range_split_at<DIM>(range); 
        auto f_regdim = function_split_and_integrate_at<DIM>(f,monte_carlo(rng,1),range_rest);

        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(bin_resolution[i]);
        std::size_t factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) factor*=bin_resolution[i];
        MutexedTensorVector<const Reg*,DIMBINS> perbin(bin_resolution,nmutexes);
        std::size_t progress = 0; std::size_t final_progress = seq_regions.size();
        auto logger_bins = logger_step(logger, "region bin stratification");
        logger_bins.log_progress(progress,final_progress);
        std::for_each(std::execution::seq, seq_regions.begin(),seq_regions.end(),[&] (const auto& r) {
            for (auto pos : pixels_in_region(r,bin_resolution,range_first)) {
                perbin.push_at(pos,&r);
            } 
        });  
        logger_bins.log_progress(final_progress,final_progress);
        auto logger_control_variates = logger_step(logger,"residual and variance reduction");
        for_each(sequential,multidimensional_range(bin_resolution),
            [&] (const std::array<std::size_t, DIMBINS>& pos) {
                using Sample = std::decay_t<decltype(f_regdim(std::declval<std::array<Float,DIM>>()))>;
                Range<Float,DIM> binrange = range_first;
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
                    const auto& [sample, sfactor] = region_sampler.sample(r,region_bin_range,rng); 
                    accumulator.push(
                        f_regdim(sample)*double(factor)*rrfactor*sfactor,
                        r->approximation_at(sample)*double(factor)*rrfactor*sfactor
                    );
//                    std::cerr<<"Usage : "<<f_regdim(sample)<<" "<<r->approximation_at(sample)<<std::endl;
                }
                bins(pos) = accumulator.integral(approximation);
//                std::cerr<<"Integral : "<<accumulator.integral(approximation);    
        }   ,logger_control_variates);
    }
};

template<typename RR, typename CV, typename RS, typename RNG>
auto regions_integrator_variance_reduction(RR&& rr, CV&& cv, RS&& rs, RNG&& rng, unsigned long samples, std::size_t nmutexes = 16, std::enable_if_t<!std::is_integral_v<RNG>,int> dummy = 0) {
    return RegionsIntegratorVarianceReduction<RR,CV,RS,RNG>(std::forward<RR>(rr), std::forward<CV>(cv), std::forward<RS>(rs), std::forward<RNG>(rng),samples,nmutexes);
}

template<typename RR, typename CV, typename RNG>
auto regions_integrator_variance_reduction(RR&& rr, CV&& cv, RNG&& rng, unsigned long samples, std::size_t nmutexes = 16, std::enable_if_t<!std::is_integral_v<RNG>,int> dummy = 0) {
    return regions_integrator_variance_reduction(std::forward<RR>(rr), std::forward<CV>(cv), region_sampling_uniform(), std::forward<RNG>(rng),samples,nmutexes);
}

template<typename RR, typename CV>
auto regions_integrator_variance_reduction(RR&& rr, CV&& cv, unsigned long samples, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return regions_integrator_variance_reduction(std::forward<RR>(rr), std::forward<CV>(cv), std::mt19937(seed),samples,nmutexes);
}

template<typename RR>
auto regions_integrator_variance_reduction(RR&& rr, double alpha, unsigned long samples, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return regions_integrator_variance_reduction(std::forward<RR>(rr), cv_fixed_weight(alpha), std::mt19937(seed),samples,nmutexes);
}

}