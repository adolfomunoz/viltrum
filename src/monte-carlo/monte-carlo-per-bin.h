#pragma once
#include <array>
#include <random>
#include "../multidimensional-range.h"
#include "../range.h"
#include "../range-infinite.h"
#include "random-sequence.h"


#if (__cplusplus < 201703L)
namespace std {
    template< class T >
    inline constexpr bool is_integral_v = std::is_integral<T>::value;
}
#endif


namespace viltrum {

template<typename RNG>
class MonteCarloPerBin {
    mutable RNG rng;
    unsigned long samples;

public:
    MonteCarloPerBin(RNG&& r, unsigned long s) : rng(std::move(r)), samples(s) {}
    //If reference (and not moved) we get a random seed
    MonteCarloPerBin(RNG& r, unsigned long s) : rng(std::size_t(r())), samples(s) {} 

    //If copy constructor, we regenerate a random seed, but this is predictable in 
    // parallel cases.
    MonteCarloPerBin(const MonteCarloPerBin& mc) :
        rng(std::size_t(mc.rng())), samples(mc.samples) {} 

    void seed(std::size_t s) { rng.seed(s); }
    RNG& random_number_generator() { return rng; }
    const RNG& random_number_generator() const { return rng; }

    //This should be more efficient than integrator_per_bin(monte_carlo(...))
	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {
        std::size_t resolution_factor(1);   
        for (std::size_t i = 0; i < DIMBINS; ++i) resolution_factor*=bin_resolution[i];
        double factor = range.volume()/double(samples);
        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/((Float)(bin_resolution[i]));
        std::size_t r = 0;

        for (auto pos : multidimensional_range(bin_resolution)) {
            logger.log_progress(r++,resolution_factor);
            Range<Float,DIM> subrange = range;
            for (std::size_t i=0;i<DIMBINS;++i)
                subrange = subrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);
            for (unsigned long s = 0; s<samples;++s) {
                std::array<Float,DIM> sample;
                for (std::size_t i=0;i<DIM;++i) {
                    std::uniform_real_distribution<Float> dis(subrange.min(i),subrange.max(i));
                    sample[i] = dis(rng);
                }
                bins(pos) += f(sample)*factor;
            }
        }
        logger.log_progress(resolution_factor,resolution_factor);
	}

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const RangeInfinite<Float>& range, Logger& logger) const {
        std::size_t resolution_factor(1);   
        for (std::size_t i = 0; i < DIMBINS; ++i) resolution_factor*=bin_resolution[i];
        double factor = range.volume()/double(samples);
        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/((Float)(bin_resolution[i]));
        std::size_t r = 0;

        for (auto pos : multidimensional_range(bin_resolution)) {
            logger.log_progress(r++,resolution_factor);
            RangeInfinite<Float> subrange = range;
            for (std::size_t i=0;i<DIMBINS;++i)
                subrange = subrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);
            for (unsigned long s = 0; s<samples;++s) {
                bins(pos) += f(random_sequence<Float,RNG>(subrange,unsigned(rng())))*factor;
                // --> Unsigned rng seeds the rng. We need to check if there is a better way
            }
        }
        logger.log_progress(resolution_factor,resolution_factor);
	}

};

template<typename RNG>
auto monte_carlo_per_bin(RNG&& rng, unsigned long samples, std::enable_if_t<!std::is_integral_v<RNG>,int> dummy = 0) {
    return MonteCarloPerBin<std::decay_t<RNG>>(std::forward<RNG>(rng),samples);
}

auto monte_carlo_per_bin(unsigned long samples, std::size_t seed = std::random_device()()) {
    return monte_carlo_per_bin(std::mt19937(seed),samples);
}

}


