#pragma once

#include <array>
#include <random>
#include "../quadrature/vector-dimensions.h"
#include "../quadrature/range.h"


class RayMarching {
    unsigned long m_nb_steps;

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
        
        Float delta_t = (range.max(0)-range.min(0))/(Float)m_nb_steps;

        Float Tr = 1;
        auto t = range.min(0)+delta_t/2.;

        while(true)
        {

            if( t>range.max(0) )
                break;

            Tr *= exp(-f(t)*delta_t); 

            
            if( verbose )
                printf("Step: %f [%f, %f] - mu_t(t)= %f : T(t) = %0.10f \n", t, range.min(0), range.max(0), f(t), Tr);

            t += delta_t;
            ++ samples.counter_queries;
        }


        //printf("End--------\n");

        samples.sumatory = Tr;
        samples.counter = 1;;
    }

    template<typename F, typename Result>
    Result integral(const F& f, const Samples<Result>& samples) const {
        return (samples.counter==0)?decltype(samples.sumatory)(0):(samples.sumatory/samples.counter);
    }

    template<typename F, typename Result>
    unsigned long queries(const F& f, const Samples<Result>& samples) const {
        return (samples.counter==0)?decltype(samples.counter_queries)(0):(samples.counter_queries);
    }

    RayMarching(unsigned long nb_steps) :
        m_nb_steps(nb_steps){ }
};


auto ray_marching(unsigned long samples) {
    return RayMarching(samples);
}

template<typename RNG>
auto ray_marching(RNG&& rng, unsigned long samples) {
    return ray_marching(samples);
}


auto ray_marching(unsigned long samples, std::size_t seed) {
    return ray_marching(samples);
}

