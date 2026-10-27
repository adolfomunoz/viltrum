#pragma once
#include <vector>
#include <array>
#include <algorithm>
#include "../range.h"
#include "mutexed-bins-aligned.h"

namespace viltrum {



/**
 * @brief This integrator is not intended to be directly used by the user. Is a particular piece
 * common to all integrators that are based on Newton-Cotes regions. 
 */
class RegionsIntegratorParallelRegions {
    std::size_t nmutexes;
public:
    RegionsIntegratorParallelRegions(std::size_t n = 16) : nmutexes(n) {}
	template<typename Bins, std::size_t DIMBINS, typename SeqRegions, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate_regions(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const SeqRegions& seq_regions, const F& f, const Range<Float,DIM>& range, Logger& logger) const {
        
        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(bin_resolution[i]);
        std::size_t factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) factor*=bin_resolution[i];
        //std::atomic<std::size_t> progress = 0; 
        std::size_t final_progress = seq_regions.size();
        MutexedBinsAligned<Bins,DIMBINS> mutexedbins(bins,bin_resolution,nmutexes);
/*        std::thread for_log([&logger,&progress,final_progress] () {
             while (std::size_t(progress)<final_progress) {
                logger.log_progress(std::size_t(progress),final_progress);
                std::this_thread::sleep_for(std::chrono::milliseconds(250));
            } 
        });*/
        logger.log_progress(std::size_t(0),final_progress);
        std::for_each(std::execution::par, seq_regions.begin(),seq_regions.end(),[&] (const auto& r) {
//            progress++;
            for (auto pos: pixels_in_region(r,bin_resolution,range)) {
                Range<Float,DIM> binrange = range;
                for (std::size_t i=0;i<DIMBINS;++i)
                    binrange = binrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);

                Range<Float,DIM> region_bin_range = binrange.intersection(r.range());
                if (!region_bin_range.empty()) mutexedbins.add(pos,double(factor)*r.integral_subrange(region_bin_range));
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