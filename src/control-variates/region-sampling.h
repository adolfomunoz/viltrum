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
        std::uniform_real_distribution<Float> dis(0,1);
        for (std::size_t i=0;i<DIM;++i) {
            sample[i] = dis(rng);
        }
        
        auto pos = reg->sample_subrange(sample,range,norm);                       
        Float pdf = reg->pdf_subrange(pos,range,norm);
        return std::tuple<std::array<Float,DIM>,Float>(pos,(pdf<1.e-10)?0.0:(1.0/pdf));
    } 
};

template<typename Norm = NormDefault>
class region_sampling_russian_roulette {
    Norm norm;
public:
    region_sampling_russian_roulette(const Norm& n = NormDefault()) : norm(n) {}
    template<typename R, typename Float, std::size_t DIM, typename RNG>
    std::tuple<std::array<Float,DIM>,Float> sample(const R& reg, const Range<Float,DIM>& range, RNG& rng) const {
        std::array<Float,DIM> sample;
        std::uniform_real_distribution<Float> dis(0,1);
        for (std::size_t i=0;i<DIM;++i) {
            sample[i] = dis(rng);
        }
        //importances for russian roulette. These probabilities do
        // not work propperly, we need to find new ones
        Float pdf_importance = reg->pdf_integral_subrange(range,norm);
        Float pdf_uniform = range.volume();
        Float rr_sample = dis(rng);

        Float pdf = 0;
        std::array<Float,DIM> pos;
        // A bit weird, we are doing a RR between sampling strategies according
        // to their perceived importance.
        if (rr_sample < (pdf_importance/(pdf_importance+pdf_uniform))) {
            // importance sampling
            pos = reg->sample_subrange(sample,range,norm);                       
            pdf = reg->pdf_subrange(pos,range,norm);
        } else {
            //uniform sampling
            for (std::size_t i=0;i<DIM;++i) {
                pos[i] = range.min(i) + sample[i]*(range.max(i)-range.min(i));
            }
            pdf = Float(1)/range.volume();
        }

        
        return std::tuple<std::array<Float,DIM>,Float>(pos,(pdf<1.e-10)?0.0:(1.0/pdf));
    } 
};

template<typename Norm = NormDefault>
class region_sampling_mis {
    double power, cutoff;
    Norm norm;
public:
    region_sampling_mis(double power = 1, double cutoff = 0, const Norm& n = NormDefault()) 
        : power(power), cutoff(cutoff), norm(n) {}
    template<typename R, typename Float, std::size_t DIM, typename RNG>
    std::tuple<std::array<Float,DIM>,Float> sample(const R& reg, const Range<Float,DIM>& range, RNG& rng) const {
        std::array<Float,DIM> sample;
        std::uniform_real_distribution<Float> dis(0,1);
        for (std::size_t i=0;i<DIM;++i) {
            sample[i] = dis(rng);
        }
            //importance sampling
        auto pos = reg->sample_subrange(sample,range,norm);                       
        Float pdf_importance = reg->pdf_subrange(pos,range,norm);
        Float pdf_uniform = Float(1)/range.volume();

        Float mis_importance = std::pow(pdf_importance, power);
        Float mis_uniform = std::pow(pdf_uniform, power);
        Float mis_sum = mis_importance + mis_uniform;

        if (mis_importance < (cutoff*mis_sum)) {
            mis_importance = 0; mis_sum = mis_uniform;
        }
        if (mis_uniform < (cutoff*mis_sum)) {
            mis_uniform = 0; mis_sum = mis_importance;
        }

        Float p_rr_importance = 
            (pdf_importance > 1.e-10)?(mis_importance/pdf_importance):Float(0);
        Float p_rr_uniform = mis_uniform*range.volume(); //pdf = 1/range.volume();
        Float p_rr_sum = p_rr_importance + p_rr_uniform;
        //With these rr probabilities, weight is actually constant
        // and unless the range volume is zero it will be non-zero
        Float weight = p_rr_sum/mis_sum;

        Float rr_sample = dis(rng);

        // A bit weird, we are doing mis and RR to select which part
        if (rr_sample < (p_rr_uniform/p_rr_sum)) {
            //uniform sampling
            for (std::size_t i=0;i<DIM;++i) {
                pos[i] = range.min(i) + sample[i]*(range.max(i)-range.min(i));
             }
        } //Otherwise we leave it as is, importance sampling

        return std::tuple<std::array<Float,DIM>,Float>(pos,weight);
    } 
};

template<typename Norm = NormDefault>
class region_sampling_mis_unstable {
    double power, cutoff;
    Norm norm;
public:
    region_sampling_mis_unstable(double power = 1, double cutoff = 0, const Norm& n = NormDefault()) : 
        power(power), cutoff(cutoff), norm(n) {}
    template<typename R, typename Float, std::size_t DIM, typename RNG>
    std::tuple<std::array<Float,DIM>,Float> sample(const R& reg, const Range<Float,DIM>& range, RNG& rng) const {
        std::array<Float,DIM> sample;
        std::uniform_real_distribution<Float> dis(0,1);
        for (std::size_t i=0;i<DIM;++i) {
            sample[i] = dis(rng);
        }
            //importance sampling
        auto pos = reg->sample_subrange(sample,range,norm);                       
        Float pdf_importance = reg->pdf_subrange(pos,range,norm);
        Float pdf_uniform = Float(1)/range.volume();

        Float mis_importance = std::pow(pdf_importance, power);
        Float mis_uniform = std::pow(pdf_uniform, power);
        Float mis_sum = mis_importance + mis_uniform;

        if (mis_importance < (cutoff*mis_sum)) {
            mis_importance = 0; mis_sum = mis_uniform;
        }
        if (mis_uniform < (cutoff*mis_sum)) {
            mis_uniform = 0; mis_sum = mis_importance;
        }

        Float rr_sample = dis(rng);

        Float weight = 0;
        // A bit weird, we are doing mis and RR to select which part
        if (rr_sample < (mis_importance/mis_sum)) {
            if (pdf_importance > 1.e-10) weight = Float(1)/pdf_importance;
        } else {
            //uniform sampling
            for (std::size_t i=0;i<DIM;++i) {
                pos[i] = range.min(i) + sample[i]*(range.max(i)-range.min(i));
             }
            weight = range.volme();
        }

        
        return std::tuple<std::array<Float,DIM>,Float>(pos,weight);
    } 
};


} 