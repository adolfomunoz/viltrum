#pragma once
#include <random>
#include <execution>
#include "../range.h"
#include "mutexed-tensor-vector.h"
#include "../foreach.h"



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

class cv_fixed_weight {
    double alpha; //This should go from zero (full importance sampling) to 1 (full control variates)
public:
    cv_fixed_weight(double a = 1) : alpha(a) {}

    template<typename Sample>
    class Accumulator {
    private:
        Sample sum; std::size_t size; double alpha;
        Accumulator(double alpha, const Sample& ini = Sample(0)) : 
            sum(ini), size(0), alpha(alpha) {}
    public:
        void push(const Sample& function_sample, const Sample& approximation_sample) {
            sum += function_sample - alpha*approximation_sample;
            ++size;
        }

        Sample integral(const Sample& approximation) const {
            return sum/double(size) + alpha*approximation;
        }
        friend class cv_fixed_weight;
    };

    template<typename Sample>
    Accumulator<Sample> accumulator(const Sample& ini = Sample(0)) const {
        return Accumulator(alpha,ini);
    }
};

template<typename Norm = NormDefault>
class cv_optimize_weight {
    Norm norm;
public:
    cv_optimize_weight(const Norm& n = Norm()) : norm(n) {}
    template<typename Sample>
    class Accumulator {
    private:
        Norm norm;
        Sample sum_f, sum_app; std::size_t size;
        //These are for online calculation of variance and covariance
        double k_f, k_app, e_f, e_ap, e_ap2, e_fap;
        Accumulator(const Norm& n, const Sample& ini = Sample(0)) : 
            norm(n),sum_f(ini),sum_app(ini),size(0),
            k_f(0), k_app(0), e_f(0), e_ap(0), e_ap2(0), e_fap(0) {}
    public:
        void push(const Sample& function_sample, const Sample& approximation_sample) {
            //For online variance and covariance (and alpha)
            if (size==0) {
                k_f = norm(function_sample); k_app = norm(approximation_sample);
            }

            e_f += (norm(function_sample) - k_f);
            e_ap += (norm(approximation_sample) - k_app);
            e_ap2 += (norm(approximation_sample) - k_app)*(norm(approximation_sample) - k_app);
            e_fap += (norm(function_sample) - k_f)*(norm(approximation_sample) - k_app);

            sum_f += function_sample; sum_app += approximation_sample; 
            ++size;
        }

        double covariance() const {
            return (e_fap - (e_f*e_ap)/double(size))/double(size-1);
        } 

        double variance() const {
            return (e_ap2 - (e_ap*e_ap)/double(size))/double(size-1);
        } 

        double alpha() const {
            auto c = std::max(0.0,covariance());
            auto v = std::max(c,variance());
            return c/v;
        } 

        Sample integral(const Sample& approximation) const {
            auto a = alpha();
            return (sum_f - a*sum_app)/double(size) + a*approximation;
        }
        friend class cv_optimize_weight;
    };

    template<typename Sample>
    Accumulator<Sample> accumulator(const Sample& ini = Sample(0)) const {
        return Accumulator<Sample>(norm,ini);
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