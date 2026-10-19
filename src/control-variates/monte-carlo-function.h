#pragma once

#include "../range.h"
#include "../range-infinite.h"
#include "../array.h"
#include <random>

namespace viltrum {
 
template<typename F, typename RNG, typename Range>
class MonteCarloFunction {
    const F& f;
    RNG& rng;
    const Range& range_mc;
public:
    MonteCarloFunction(const F& f, RNG& rng, const Range& range_mc) :
        f(f), rng(rng), range_mc(range_mc) {}

    template<typename Float, std::size_t DIM>
    auto operator()(const std::array<Float,DIM>& x) const {
        std::array<Float, Range::DIM> sample;
        for (std::size_t i=0;i<Range::DIM;++i) {
		    std::uniform_real_distribution<Float> dis(range_mc.min(i),range_mc.max(i));
		    sample[i] = dis(rng);
	    }
        return f(x | sample);
    }
};

template<typename F, typename RNG, typename Float>
class MonteCarloFunction<F,RNG,RangeInfinite<Float>> {
    const F& f;
    RNG& rng;
    const RangeInfinite<Float>& range_mc;
public:
    MonteCarloFunction(const F& f, RNG& rng, const RangeInfinite<Float>& range_mc) :
        f(f), rng(rng), range_mc(range_mc) {}
    template<std::size_t DIM>
    auto operator()(const std::array<Float,DIM>& x) const {
        return f(concat(x,random_sequence(range_mc,rng)));
    }
};

template<typename F, typename RNG, typename Range>
auto monte_carlo_function(const F& f, RNG& rng, const Range& range) {
    return MonteCarloFunction<F,RNG,Range>(f,rng,range);
}

}