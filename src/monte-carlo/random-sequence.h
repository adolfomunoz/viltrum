#pragma once
#include <random>
#include "random-sequence-seed.h"
#include "random-sequence-rng.h"
#include "random-sequence-ref.h"




namespace viltrum {

template<typename Number,typename RNG = std::mt19937>
auto random_sequence(const RangeInfinite<Number>& r, unsigned seed = std::random_device()()) {
    return random_sequence_rng(r,seed);
} 

template<typename Number,typename RNG>
auto random_sequence(const RangeInfinite<Number>& r, RNG& rng) {
    return random_sequence_ref(r,rng);
} 



}
