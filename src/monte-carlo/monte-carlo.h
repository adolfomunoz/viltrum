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
class MonteCarlo {
    mutable RNG rng;
    unsigned long samples;

public:
    MonteCarlo(RNG&& r, unsigned long s) : rng(std::move(r)), samples(s) {}
    //If reference (and not moved) we get a random seed
    MonteCarlo(RNG& r, unsigned long s) : rng(std::size_t(r())), samples(s) {} 

    //If copy constructor, we regenerate a random seed, but this is predictable in 
    // parallel cases.
    MonteCarlo(const MonteCarlo& mc) :
        rng(std::size_t(mc.rng())), samples(mc.samples) {} 

    void seed(std::size_t s) { rng.seed(s); }
    RNG& random_number_generator() { return rng; }
    const RNG& random_number_generator() const { return rng; }

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {

        double resolution_factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) resolution_factor*=bin_resolution[i];
        double factor = resolution_factor*range.volume()/double(samples);
        unsigned long i;
        for (i=0;i<samples;++i) {
            logger.log_progress(i,samples);
            std::array<Float,DIM> sample;
	        for (std::size_t i=0;i<DIM;++i) {
		        std::uniform_real_distribution<Float> dis(range.min(i),range.max(i));
		        sample[i] = dis(rng);
	        }
            if (range.is_inside(sample)) {
                std::array<std::size_t,DIMBINS> pos;
                for (std::size_t i=0;i<DIMBINS;++i) {
                    pos[i] = std::size_t(bin_resolution[i]*(sample[i] - range.min(i))/(range.max(i) - range.min(i)));
                }
                bins(pos) += f(sample)*factor;
            }
        }
        logger.log_progress(samples,samples);
	}

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const RangeInfinite<Float>& range, Logger& logger) const {

        double resolution_factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) resolution_factor*=bin_resolution[i];
        double factor = resolution_factor*range.volume()/double(samples);
        unsigned long i;
        for (i=0;i<samples;++i) {
            logger.log_progress(i,samples);
            auto sample = random_sequence<Float,RNG>(range,std::size_t(rng())); //Seeds the RNG
            std::array<std::size_t,DIMBINS> pos;
            auto it = sample.begin();
            for (std::size_t i=0;i<DIMBINS;++i,++it) {
                pos[i] = std::size_t(bin_resolution[i]*((*it) - range.min(i))/(range.max(i) - range.min(i)));
            }
            bins(pos) += f(sample)*factor;
        }
        logger.log_progress(samples,samples);
	}

};

template<typename RNG>
auto monte_carlo(RNG&& rng, unsigned long samples, std::enable_if_t<!std::is_integral_v<RNG>,int> dummy = 0) {
    return MonteCarlo<std::decay_t<RNG>>(std::forward<RNG>(rng),samples);
}

auto monte_carlo(unsigned long samples, std::size_t seed = std::random_device()()) {
    return monte_carlo(std::mt19937(seed),samples);
}

}


