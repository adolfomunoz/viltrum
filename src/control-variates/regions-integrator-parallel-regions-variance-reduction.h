#pragma once
#include <vector>
#include <cstdint>
#include <cassert>
#include <mutex>
#include <atomic>
#include <thread>
#include <chrono>
#include "../range.h"
#include "../newton-cotes/mutexed-bins-aligned.h"
#include "../newton-cotes/region-bin-accumulator.h"
#include "../foreach.h"
#include "../progress-logger.h"
#include "region-stratification.h"
#include "region-sampling.h"

namespace viltrum {


/**
 * @brief This integrator is not intended to be directly used by the user. It is a variance reduction
 * region integrator distributes the samples and parallelizes among all regions of the whole
 * integrand, independently of the number and distribution of bins.
 */
template<typename Stratification, typename RS, typename RNG>
class RegionsIntegratorParallelRegionsVarianceReduction {
    Stratification stratification;
    RS region_sampler;
    mutable RNG rng;
    unsigned long samples;
    double alpha; //Weight between CV and residual. 0 means only CV, 1 means only residual.
    std::size_t nmutexes;

    // SplitMix64 finalizer applied to (seed, index): a cheap, well-mixed, deterministic per-region seed.
    static std::uint64_t hash_seed(std::uint64_t seed, std::uint64_t index) {
        std::uint64_t z = seed + (index + 1) * 0x9E3779B97F4A7C15ull;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
public:

    RegionsIntegratorParallelRegionsVarianceReduction(Stratification&& stratification, RS&& rs, RNG&& r, unsigned long s, double a = 1.0, std::size_t n = 16) 
        : stratification(stratification), region_sampler(rs), rng(std::move(r)), samples(s), alpha(a), nmutexes(n) {}
    RegionsIntegratorParallelRegionsVarianceReduction(Stratification&& stratification, RS&& rs, RNG& r, unsigned long s, double a = 1.0, std::size_t n = 16) 
        : stratification(stratification), region_sampler(rs), rng(std::size_t(r())), samples(s), alpha(a), nmutexes(n) {}	
    
    template<typename Bins, std::size_t DIMBINS, typename SeqRegions, typename F, typename IntegrationRange, typename Logger>
	void integrate_regions(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const SeqRegions& seq_regions, const F& f, const IntegrationRange& range, Logger& logger) const {
        
        using Reg = typename SeqRegions::value_type;
        using RegRange = std::decay_t<decltype(std::declval<Reg>().range())>;
        constexpr std::size_t DIM = RegRange::size;
        using Float = typename RegRange::value_type;
        auto ranges = range_split_at<DIM>(range);
        auto range_first = std::get<0>(ranges); auto range_rest = std::get<1>(ranges);

        //Random number generators for each region, so the seed makes this predictable.


        //Instead of creating a new random number per region, we hash a seed for each
        // region lazily, when needed. This avods so many seeds and the associated memory
        // storage.
        std::size_t rng_seed = std::size_t(rng());

        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(bin_resolution[i]);
        std::size_t factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) factor*=bin_resolution[i];
        unsigned long total_samples = samples*factor; //Total samples is the number of samples per bin times the number of bins
        //We should double check that there is no overflow in the multiplication above, but it should be ok for now.
        
        auto samples_per_region = stratification.samples_per_region(total_samples, seq_regions, rng);

        //We index the regions so we can access all their information in parallel
        std::vector<std::size_t> region_index(seq_regions.size());
        std::iota(region_index.begin(),region_index.end(),0);

        std::atomic<std::size_t> progress{0};
        std::size_t final_progress = seq_regions.size();
        MutexedBinsAligned<Bins,DIMBINS> mutexedbins(bins,bin_resolution,nmutexes);

        // Only the helper thread calls the logger; workers just bump a relaxed counter.
        auto for_log = progress_logger(logger,final_progress,
            [&progress] { return progress.load(std::memory_order_relaxed); });
        //We parallelize the regions using "par" because "par_unseq" does not allow for mutexes inside.
        std::for_each(std::execution::par, region_index.begin(),region_index.end(),
            [&] (const std::size_t& i) {
                const auto& r = seq_regions[i];
                unsigned long samples = samples_per_region[i];

                //The contributions of this region (the control variate integral in each bin and one for each
                // sample of the residual) are accumulated locally and added to the shared bins once, at the end. 
                // Adding each sample directly takes a lock per sample on the few bins that the region overlaps, 
                // and with samples this cheap the threads spend most of the time waiting for those locks.
                auto bin_box = pixels_in_region(r,bin_resolution,range_first);
                // The type of the contributions, not the one returned by bins(), which can be a proxy (e.g. a pixel of an image) 
                using BinValue = std::decay_t<decltype(double(factor)*alpha*r.integral_subrange(r.range()))>;
                RegionBinAccumulator<MutexedBinsAligned<Bins,DIMBINS>,DIMBINS,BinValue> accumulator(mutexedbins,bin_box);

                for (auto pos: bin_box) {
                    Range<Float,DIM> binrange = r.range();
                    for (std::size_t i=0;i<DIMBINS;++i)
                        binrange = binrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);

                    Range<Float,DIM> region_bin_range = binrange.intersection(r.range());
                    if (!region_bin_range.empty()) accumulator.add(pos,double(factor)*alpha*r.integral_subrange(region_bin_range));
                }
                
                //Now we sample the residual, and distribute the samples according to the stratification strategy
                if (samples > 0) { 
                    // The RNG must be declared before f_regdim, which keeps a reference to it
                    RNG rng(static_cast<typename RNG::result_type>(hash_seed(rng_seed,i)));
                    auto f_regdim = function_split_and_integrate_at<DIM>(f,monte_carlo_per_bin(rng,1),range_rest);
                    double region_factor = double(factor)/double(samples);

                    for (unsigned long s = 0; s<samples; ++s) {
                        const auto& [sample, sfactor] = region_sampler.sample(&r,r.range(),rng);
                        assert(r.range().is_inside(sample));
                        std::array<std::size_t,DIMBINS> pos;
                        for (std::size_t i=0;i<DIMBINS;++i) {
                            pos[i] = std::size_t(bin_resolution[i]*(sample[i] - range.min(i))/(range.max(i) - range.min(i)));
                            if (pos[i]>=bin_resolution[i]) { //This happens very few times, when the sample gets to exactly 1
                                pos[i] = bin_resolution[i]-1; 
                            }
                        }

                        accumulator.add(pos,
                            region_factor*sfactor*
                                (f_regdim(sample) - alpha*r.approximation_at(sample))
                            );
                    } 
                }
                accumulator.flush();
                progress.fetch_add(1,std::memory_order_relaxed);
            }
        );
        
        for_log.finish();
    }
};


/**
 * @brief This integrator is not intended to be directly used by the user. Is a particular piece
 * common to all integrators that are based on Newton-Cotes regions. The integrand is therefore ignored
 * 
 * @tparam Stratification The type of the stratification strategy.
 * @tparam RS The type of the region sampler.
 * @tparam RNG The type of the random number generator.
 * @tparam std::size_t Number of mutexes. The more (allowed by the hardware and OS) the less conflicts.
 */
template<typename Stratification, typename RS, typename RNG>
auto regions_integrator_parallel_regions_variance_reduction(
    Stratification&& stratification, RS&& rs, RNG&& rng, unsigned long samples_per_bin, double alpha = 1.0, std::size_t nmutexes = 1024, std::enable_if_t<!std::is_integral_v<RNG>,int> dummy = 0) {
    return RegionsIntegratorParallelRegionsVarianceReduction<
        std::decay_t<Stratification>, std::decay_t<RS>, std::decay_t<RNG>>(
            std::forward<Stratification>(stratification),
            std::forward<RS>(rs),
            std::forward<RNG>(rng),
            samples_per_bin,
            alpha,
            nmutexes);
}

template<typename Stratification, typename RS>
auto regions_integrator_parallel_regions_variance_reduction(
    Stratification&& stratification, RS&& rs, unsigned long samples_per_bin, std::size_t seed = std::random_device()(), double alpha = 1.0, std::size_t nmutexes = 1024) {
    return regions_integrator_parallel_regions_variance_reduction(
        std::forward<Stratification>(stratification),
        std::forward<RS>(rs),
        std::mt19937(seed),
        samples_per_bin,
        alpha,
        nmutexes);
    }


}