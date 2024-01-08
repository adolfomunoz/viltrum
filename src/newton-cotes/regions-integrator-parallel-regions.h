#pragma once
#include <vector>
#include <mutex>
#include "../range.h"

namespace viltrum {

template<typename Bins, std::size_t DIMBINS>
class MutexedBins {
    Bins& bins;
    std::array<std::size_t,DIMBINS> bin_resolution;
    std::vector<std::mutex> mutexes;
    using value_type = std::decay_t<decltype(bins(std::array<std::size_t,DIMBINS>{}))>;
    constexpr std::size_t hash_of(const std::array<std::size_t,DIMBINS>& index) {
        std::size_t sum{0}, prod{1};
        for (std::size_t i = 0; i<DIMBINS; ++i) { sum+=prod*index[i]; prod*=bin_resolution[i]; }
        return sum % mutexes.size();
    }
public:
    MutexedBins(Bins& b, const std::array<std::size_t,DIMBINS>& br,std::size_t nmutexes) : 
        bins(b),bin_resolution(br),mutexes(nmutexes) {}
    void add(const std::array<std::size_t,DIMBINS>& i,const value_type& v) {
        std::scoped_lock lock(mutexes[hash_of(i)]);
        bins(i) += v;
    }
};

/**
 * @brief This integrator is not intended to be directly used by the user. Is a particular piece
 * common to all integrators that are based on Newton-Cotes regions. 
 */
class RegionsIntegratorParallelRegions {
    std::size_t nmutexes;
public:
    RegionsIntegratorParallelRegions(std::size_t n = 16) : nmutexes(n) {}
	template<typename Bins, std::size_t DIMBINS, typename SeqRegions, typename Float, std::size_t DIM, typename Logger>
	void integrate_regions(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const SeqRegions& seq_regions, const Range<Float,DIM>& range, Logger& logger) const {
        
        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(bin_resolution[i]);
        std::size_t factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) factor*=bin_resolution[i];
        std::atomic<std::size_t> progress = 0; std::size_t final_progress = seq_regions.size();
        MutexedBins<Bins,DIMBINS> mutexedbins(bins,bin_resolution,nmutexes);
        std::thread for_log([&logger,&progress,final_progress] () {
             while (std::size_t(progress)<final_progress) {
                logger.log_progress(std::size_t(progress),final_progress);
                std::this_thread::sleep_for(std::chrono::milliseconds(250));
            } 
        });
        std::for_each(std::execution::par_unseq, seq_regions.begin(),seq_regions.end(),[&] (const auto& r) {
            progress++;
            std::array<std::size_t,DIMBINS> start_bin, end_bin;
            for (std::size_t i = 0; i<DIMBINS;++i) {
                start_bin[i] = std::size_t(std::max(Float(0),(r.range().min(i) - range.min(i))/drange[i]));
                end_bin[i]   = std::min(bin_resolution[i],std::size_t((r.range().max(i) - range.min(i))/drange[i])+1);
            }
            if (start_bin == end_bin) mutexedbins.add(start_bin,r.integral());
            else for (auto pos : multidimensional_range(start_bin, end_bin)) {
                Range<Float,DIM> binrange = range;
                for (std::size_t i=0;i<DIMBINS;++i)
                    binrange = binrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);

                Range<Float,DIM> region_bin_range = binrange.intersection(r.range());
                if (!region_bin_range.empty()) mutexedbins.add(pos,double(factor)*r.integral_subrange(region_bin_range));
            } 
        });
        
        for_log.join();  
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