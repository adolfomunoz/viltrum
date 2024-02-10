#pragma once
#include <random>
#include "../range.h"
#include "mutexed-tensor-vector.h"

namespace viltrum {

class RRUniformRegion {
    //This is the MonteCarlo residual, RR among regions
    std::uniform_int_distribution<std::size_t> rr;
public:
    template<typename Rs>
    RRUniformRegion(const Rs& regions) : rr(0,regions.size()-1) {}

    template<typename RNG>
    std::tuple<std::size_t,double> choose(RNG& rng) { 
        return std::tuple(rr(rng),rr.max()+1); 
    }
};

class RRIntegralRegion {
    //This is the MonteCarlo residual, RR among regions
    std::discrete_distribution<std::size_t> rr;
    std::vector<double> weights;
    float norm(float v) const { return v; }
    double norm(double v) const { return v; }
    template<typename V>
    auto norm(const V& v, typename std::enable_if_t<std::is_arithmetic_v<typename V::value_type>, int> = 0) const {
        auto i = v.begin(); 
        using S = decltype(norm(*i));
        S s; bool first = true;
        if (i != v.end()) { s = norm(*i); ++i; }
        while (i != v.end()) { s += norm(*i); ++i; }
        return s;
    }     
    
public:
    template<typename Rs>
    RRIntegralRegion(const Rs& regions) : weights(regions.size()) {
        std::size_t i = 0;
        for (const auto& [r,region_bin_range] : regions) {
            weights[i++] = norm(r->integral_subrange(region_bin_range));
	    }
        rr = std::discrete_distribution<std::size_t>(weights.begin(),weights.end());
    }

    template<typename RNG>
    std::tuple<std::size_t,double> choose(RNG& rng) { 
        std::size_t choice = rr(rng);
        return std::tuple(choice,1.0/rr.probabilities()[choice]); 
    }
};



template<typename RR, typename RNG>
class RegionsIntegratorParallelVarianceReduction {
    mutable RNG rng;
    unsigned long samples;
    std::size_t nmutexes;

public:
    RegionsIntegratorParallelVarianceReduction(RNG&& r, unsigned long s,std::size_t n = 16) 
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
        MutexedTensorVector<std::tuple<const Reg*,Range<Float,DIM>>,DIMBINS> perbin(bin_resolution,nmutexes);
        auto logger_bins = logger_step(logger, "region bin stratification");
        logger_bins.log_progress(progress,final_progress);
        std::for_each(std::execution::par_unseq, seq_regions.begin(),seq_regions.end(),[&] (const auto& r) {
            for (auto pos : pixels_in_region(r,bin_resolution,range)) {
                Range<Float,DIM> binrange = range;
                for (std::size_t i=0;i<DIMBINS;++i)
                    binrange = binrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);

                Range<Float,DIM> region_bin_range = binrange.intersection(r.range());
                if (!region_bin_range.empty()) perbin.push_at(pos,std::tuple(&r,region_bin_range));
            } 
        });  
        logger_bins.log_progress(final_progress,final_progress);
        auto logger_control_variates = logger_step(logger,"control variates");
        for_each(parallel,multidimensional_range(bin_resolution),
            [&] (const std::array<std::size_t, DIMBINS>& pos) {
                //This is the control variate
                for (const auto& [r,region_bin_range] : perbin[pos]) {
                    bins(pos) += double(factor)*r->integral_subrange(region_bin_range);
	            }
                //This is the MonteCarlo residual, RR among regions
                RR rr(perbin[pos]);
                for (unsigned long i=0;i<samples;++i) {
                    auto [chosen,rrfactor] = rr.choose(rng);
                    const auto& [r,region_bin_range] = perbin[pos][chosen];
                    std::array<Float,DIM> sample;
	                for (std::size_t i=0;i<DIM;++i) {
		                std::uniform_real_distribution<Float> dis(region_bin_range.min(i),region_bin_range.max(i));
		                sample[i] = dis(rng);
                    }
                    bins(pos) += (f(sample)-r->approximation_at(sample))*double(factor)*region_bin_range.volume()*rrfactor/double(samples);
                }
        }   ,logger_control_variates);
    }
};

template<typename RR, typename RNG>
auto regions_integrator_parallel_variance_reduction(RNG&& rng, unsigned long samples, std::size_t nmutexes = 16, std::enable_if_t<!std::is_integral_v<RNG>,int> dummy = 0) {
    return RegionsIntegratorParallelVarianceReduction<RR,RNG>(std::forward<RNG>(rng),samples,nmutexes);
}

template<typename RR>
auto regions_integrator_parallel_variance_reduction(unsigned long samples, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return regions_integrator_parallel_variance_reduction<RR>(std::mt19937(seed),samples,nmutexes);
}

}