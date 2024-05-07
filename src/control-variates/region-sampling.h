#pragma once
#include "../range.h"
#include "norm.h"
#include <array>
#include <random>

namespace viltrum {

class region_sampling_uniform {
public:
    template<typename R, typename Float, std::size_t DIM, typename RNG>
    std::tuple<std::array<Float,DIM>,Float> sample(const R& reg, const Range<Float,DIM>& range, RNG& rng) const {
        std::array<Float,DIM> sample;
        for (std::size_t i=0;i<DIM;++i) {
            std::uniform_real_distribution<Float> dis(range.min(i),range.max(i));
            sample[i] = dis(rng);
        }
        return std::tuple<std::array<Float,DIM>,Float>(sample,range.volume());
    } 
};

template<typename Norm = NormDefault>
class region_sampling_importance {
    Norm norm;
public:
    region_sampling_importance(const Norm& n = NormDefault()) : norm(n) {}
    template<typename R, typename Float, std::size_t DIM, typename RNG>
    std::tuple<std::array<Float,DIM>,Float> sample(const R& reg, const Range<Float,DIM>& range, RNG& rng) const {
        std::array<Float,DIM> sample;
        for (std::size_t i=0;i<DIM;++i) {
            std::uniform_real_distribution<Float> dis(range.min(i),range.max(i));
            sample[i] = dis(rng);
        }
        auto pos = reg->sample_subrange(sample,range,norm);
        Float den = reg->pdf_subrange(pos,range,norm);
        return std::tuple<std::array<Float,DIM>,Float>(pos,(den<1.e-10)?Float(0):range.volume()/den);
    } 
};


} 