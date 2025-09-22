#pragma once
#include <random>
#include "random-sequence-seed.h"
#include "random-sequence-rng.h"


namespace viltrum {

template<typename Number,typename RNG = std::mt19937>
auto random_sequence(const RangeInfinite<Number>& r, std::size_t seed = std::random_device()()) {
    return random_sequence_rng(r,seed);
} 




}
