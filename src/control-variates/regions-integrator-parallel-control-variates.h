#pragma once
#include <vector>
#include <mutex>
#include "../range.h"
#include "../tensor.h"





namespace viltrum {

template<typename T, std::size_t DIMBINS>
class MutexedTensorVector {
    tensor<std::vector<T>> data;
    std::vector<std::mutex> mutexes;
    constexpr std::size_t hash_of(const std::array<std::size_t,DIMBINS>& index) {
        std::size_t sum{0}, prod{1};
        for (std::size_t i = 0; i<DIMBINS; ++i) { sum+=prod*index[i]; prod*=data.resolution(i); }
        return sum % mutexes.size();
    }
public:
    MutexedTensorVector(const std::array<std::size_t, DIMBINS>& r,std::size_t nmutexes) : 
        data(b),mutexes(nmutexes) {}
    void push_at(const std::array<std::size_t,DIMBINS>& i,const T& v) {
        std::scoped_lock lock(mutexes[hash_of(i)]);
        data[i].push_back(v);
    }
    //Read only access
    const std::vector<T>& operator[](const std::array<std::size_t,DIMBINS>& p) const {
        return data[p];
    }

};

class RegionsIntegratorParallelControlVariates {
    unsigned long samples;
    std::size_t nmutexes;

public:
    RegionsIntegratorParallelControlVariates(unsigned long s,std::size_t n = 16) 
        : samples(s), nmutexes(n) {}
	template<typename Bins, std::size_t DIMBINS, typename SeqRegions, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate_regions(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const SeqRegions& seq_regions, const F& f, const Range<Float,DIM>& range, Logger& logger) const {
        
        using Reg = SeqRegions::value_type;

        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(bin_resolution[i]);
        std::size_t factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) factor*=bin_resolution[i];
        std::atomic<std::size_t> progress = 0; std::size_t final_progress = seq_regions.size();
        MutexedTensorVector<std::tuple<Reg,Range<Float,DIM>>,DIMBINS> perbin;
        std::for_each(std::execution::par_unseq, seq_regions.begin(),seq_regions.end(),[&] (const auto& r) {

            std::array<std::size_t,DIMBINS> start_bin, end_bin;
            for (std::size_t i = 0; i<DIMBINS;++i) {
                start_bin[i] = std::size_t(std::max(Float(0),(r.range().min(i) - range.min(i))/drange[i]));
                end_bin[i]   = std::min(bin_resolution[i],std::size_t((r.range().max(i) - range.min(i))/drange[i])+1);
            }
            if (start_bin == end_bin) perbin.push_at(startbin,{r,r.range()});
            else for (auto pos : multidimensional_range(start_bin, end_bin)) {
                Range<Float,DIM> binrange = range;
                for (std::size_t i=0;i<DIMBINS;++i)
                    binrange = binrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);

                Range<Float,DIM> region_bin_range = binrange.intersection(r.range());
                if (!region_bin_range.empty()) perbin.push_at(pos,{r,region_bin_range});
            } 
        });
        
//        for_log.join();  
        logger.log_progress(final_progress,final_progress);
    }
};


/**
 * @brief This integrator is not intended to be directly used by the user. Is a particular piece
 * common to all integrators that are based on Newton-Cotes regions. The integrand is therefore ignored
 * 
 * @tparam std::size_t Number of mutexes. The more (allowed by the hardware and OS) the less conflicts.
 */
auto regions_integrator_parallel_regions(std::size_t nmutexes = 16) {
    return RegionsIntegratorParallelRegions(nmutexes);
}

}