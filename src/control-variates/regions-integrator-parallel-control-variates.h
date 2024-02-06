#pragma once
#include <vector>
#include <mutex>
#include <random>
#include "../range.h"
#include "../tensor.h"





namespace viltrum {

template<typename T, std::size_t DIMBINS>
class MutexedTensorVector {
    tensor<std::vector<T>,DIMBINS> data;
    std::vector<std::mutex> mutexes;
    constexpr std::size_t hash_of(const std::array<std::size_t,DIMBINS>& index) {
        std::size_t sum{0}, prod{1};
        for (std::size_t i = 0; i<DIMBINS; ++i) { sum+=prod*index[i]; prod*=data.resolution(i); }
        return sum % mutexes.size();
    }
public:
    MutexedTensorVector(const std::array<std::size_t, DIMBINS>& r,std::size_t nmutexes) : 
        data(r), mutexes(nmutexes) {}

    void push_at(const std::array<std::size_t,DIMBINS>& i,T&& v) {
        std::scoped_lock lock(mutexes[hash_of(i)]);
        data[i].push_back(std::forward<T>(v));
    }
    //Read only access
    const std::vector<T>& operator[](const std::array<std::size_t,DIMBINS>& p) const {
        return data[p];
    }

};

template<typename RNG>
class RegionsIntegratorParallelControlVariates {
    mutable RNG rng;
    unsigned long samples;
    std::size_t nmutexes;

public:
    RegionsIntegratorParallelControlVariates(RNG&& r, unsigned long s,std::size_t n = 16) 
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
                std::uniform_int_distribution<std::size_t> rr(0,perbin[pos].size()-1);
                for (unsigned long i=0;i<samples;++i) {
                    std::size_t chosen = rr(rng);
                    const auto& [r,region_bin_range] = perbin[pos][chosen];
                    std::array<Float,DIM> sample;
	                for (std::size_t i=0;i<DIM;++i) {
		                std::uniform_real_distribution<Float> dis(region_bin_range.min(i),region_bin_range.max(i));
		                sample[i] = dis(rng);
                    }
                    bins(pos) += (f(sample)-r->approximation_at(sample))*double(factor)*region_bin_range.volume()*double(perbin[pos].size())/double(samples);
                }
        }   ,logger_control_variates);
    }
};

template<typename RNG>
auto regions_integrator_parallel_control_variates(RNG&& rng, unsigned long samples, std::size_t nmutexes = 16, std::enable_if_t<!std::is_integral_v<RNG>,int> dummy = 0) {
    return RegionsIntegratorParallelControlVariates<RNG>(std::forward<RNG>(rng),samples,nmutexes);
}

auto regions_integrator_parallel_control_variates(unsigned long samples, std::size_t seed = std::random_device()(), std::size_t nmutexes = 16) {
    return regions_integrator_parallel_control_variates(std::mt19937(seed),samples,nmutexes);
}

}