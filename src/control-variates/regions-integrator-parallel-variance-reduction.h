#pragma once
#include <random>
#include "../range.h"
#include "mutexed-tensor-vector.h"

namespace viltrum {

struct NormDefault {
    float operator()(float v) const { return v; }
    double operator()(double v) const { return v; }
    template<typename V>
    auto operator()(const V& v, typename std::enable_if_t<std::is_arithmetic_v<typename V::value_type>, int> = 0) const {
        auto i = v.begin(); 
        using S = decltype(norm(*i));
        S s; bool first = true;
        if (i != v.end()) { s = norm(*i); ++i; }
        while (i != v.end()) { s += norm(*i); ++i; }
        return s;
    }     
};

class rr_uniform_region {
public:
    class RR {
    private:
        //This is the MonteCarlo residual, RR among regions
        std::uniform_int_distribution<std::size_t> rr;
        template<typename Rs>
        RR(const Rs& regions) : rr(0,regions.size()-1) {}
    public:
        template<typename RNG>
        std::tuple<std::size_t,double> choose(RNG& rng) { 
            return std::tuple(rr(rng),rr.max()+1); 
        }
        friend class rr_uniform_region;
    };
    template<typename Regions>
    RR russian_roulette(const Regions& regions) const {
        return RR(regions);
    }
};

template<typename Norm = NormDefault>
class rr_integral_region {
    Norm norm;
public:
    rr_integral_region(const Norm& n = Norm()) : norm(n) {}
    class RR {
    private:
        //This is the MonteCarlo residual, RR among regions
        std::discrete_distribution<std::size_t> rr;
        std::vector<double> weights;
        Norm norm;

        template<typename Rs>
        RR(const Rs& regions, const Norm& n) : weights(regions.size()), norm(n) {
            std::size_t i = 0;
            for (const auto& [r,region_bin_range] : regions) {
                weights[i++] = norm(r->integral_subrange(region_bin_range));
            }
            rr = std::discrete_distribution<std::size_t>(weights.begin(),weights.end());
        }
        
    public:
        template<typename RNG>
        std::tuple<std::size_t,double> choose(RNG& rng) { 
            std::size_t choice = rr(rng);
            return std::tuple(choice,1.0/rr.probabilities()[choice]); 
        }
        friend class rr_integral_region;
    };
    template<typename Regions>
    RR russian_roulette(const Regions& regions) const {
        return RR(regions,norm);
    }
};

class control_variates_fixed_weight {
    double alpha; //This should go from zero (full importance sampling) to 1 (full control variates)
public:
    control_variates_fixed_weight(double a = 1) : alpha(a) {}

    template<typename Samples>
    auto integral(const Samples& function_samples, const Samples& approximation_samples, 
            const typename Samples::value_type& approximation) const {

        std::size_t size = 0;
        typename Samples::value_type sum(0);
        //We assume the same number of elements in both samples
        for (const auto& s : function_samples) { sum+=s; ++size; }
        for (const auto& s : approximation_samples) sum -= alpha*s;
        return sum/double(size) + alpha*approximation;
    }
};



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
        auto logger_control_variates = logger_step(logger,"residual and variance reduction");
        for_each(parallel,multidimensional_range(bin_resolution),
            [&] (const std::array<std::size_t, DIMBINS>& pos) {
                using Sample = decltype(f(std::declval<std::array<Float,DIM>>()));

                //This is the control variate
                Sample approximation(0);
                for (const auto& [r,region_bin_range] : perbin[pos]) {
                    approximation += double(factor)*r->integral_subrange(region_bin_range);
	            }
                //These are the samples for the residual
                std::list<Sample> samples_f, samples_cv;
                auto roulette = rr.russian_roulette(perbin[pos]);
                for (unsigned long i=0;i<samples;++i) {
                    auto [chosen,rrfactor] = roulette.choose(rng);
                    const auto& [r,region_bin_range] = perbin[pos][chosen];
                    std::array<Float,DIM> sample;
	                for (std::size_t i=0;i<DIM;++i) {
		                std::uniform_real_distribution<Float> dis(region_bin_range.min(i),region_bin_range.max(i));
		                sample[i] = dis(rng);
                    }
                    samples_f.push_back(f(sample)*double(factor)*rrfactor*region_bin_range.volume());
                    samples_cv.push_back(r->approximation_at(sample)*double(factor)*rrfactor*region_bin_range.volume());
                    //This is the MonteCarlo residual, RR among regions
                    bins(pos) = cv.integral(samples_f,samples_cv,approximation);
                }
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
    return regions_integrator_parallel_variance_reduction(std::forward<RR>(rr), control_variates_fixed_weight(alpha), std::mt19937(seed),samples,nmutexes);
}

}