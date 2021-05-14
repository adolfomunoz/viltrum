#pragma once

#include <array>
#include <random>
#include "../quadrature/vector-dimensions.h"
#include "../quadrature/range.h"

template<typename RNG>
class ResidualRatioTracking {
    mutable RNG rng;

    double maj_sigma_t;
    double c_sigma_t;

    template<typename Result>
    struct Samples {
        Result sumatory;
        unsigned long counter;
        unsigned long counter_queries;
        Samples() : sumatory(0),counter(0),counter_queries(0) { }
    };
public:
    template<typename F, typename Float>
    auto init(const F& f, const Range<Float,1>& range) const {
        //First element of tuple is sum, second element is number of samples
        return Samples<decltype(f(range.min(0)))>();
    }

    template<typename F, typename Float, typename Result>
    void step(const F& f, const Range<Float,1>& range, Samples<Result>& samples, bool verbose = false) const {

        std::exponential_distribution<Float> dis(maj_sigma_t);

        Float Tc = exp(-(range.max(0)-range.min(0))*c_sigma_t);

        Float Tr = 1;
        Float t = range.min(0);

        while(true)
        {
            Float dt = dis(rng);
            t += dt;

            if( t>range.max(0) )
                break;

            Tr *= (1 - (f(t) - c_sigma_t)/maj_sigma_t); 

            if( verbose )
                printf("Step: %f [%f, %f] - mu_t(t)= %f : T(t) = %f; Tc(t) = %f; T_c = %f\n", t, range.min(0), range.max(0), (f(t)-c_sigma_t)/maj_sigma_t, Tr, exp(-c_sigma_t*t), Tc);
//                printf("Step: %f [%f, %f] - mu_t(t)= %f : T(t) = %f \n", t, range.min(0), range.max(0), f(t)-c_sigma_t, Tr*exp(-t*c_sigma_t));

            ++ samples.counter_queries;
        }


        samples.sumatory += Tr*Tc;
        ++samples.counter;
    }

    template<typename F, typename Result>
    Result integral(const F& f, const Samples<Result>& samples) const {
        return (samples.counter==0)?decltype(samples.sumatory)(0):(samples.sumatory/samples.counter);
    }

    template<typename F, typename Result>
    unsigned long queries(const F& f, const Samples<Result>& samples) const {
        return (samples.counter==0)?decltype(samples.counter_queries)(0):(samples.counter_queries);
    }

    ResidualRatioTracking(RNG&& r, double majorant, double control) :
        rng(std::forward<RNG>(r)), maj_sigma_t(majorant), c_sigma_t(control) { }
};


template<typename RNG>
auto ratio_tracking(RNG&& rng, double majorant) {
    return ResidualRatioTracking(std::forward<RNG>(rng), majorant, 0.);
}

template<typename RNG>
auto residual_ratio_tracking(RNG&& rng, double majorant, double control) {
    return ResidualRatioTracking(std::forward<RNG>(rng), majorant, control);
}

/*auto ratio_tracking(double majorant, double control) {
    return ratio_tracking(std::mt19937_64(seed), majorant, control);
}*/
