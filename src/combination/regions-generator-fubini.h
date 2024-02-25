#pragma once

#include "fubini.h"

namespace viltrum {

template<typename BaseGenerator, typename RestIntegrator, std::size_t N>
class RegionsGeneratorFubini {
    BaseGenerator base_generator;
    RestIntegrator rest_integrator;

public:
    RegionsGeneratorFubini(const BaseGenerator& bg, const RestIntegrator& rs) 
        : base_generator(bg), rest_integrator(rs) {}

    template<std::size_t DIMBINS, typename F, typename IntegrationRange, typename Logger>
    auto generate(const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const IntegrationRange& range, Logger& logger) const {
            auto [range_first, range_rest] = range_split_at<N>(range);
            auto f_gen = function_split_and_integrate_at<N>(f,rest_integrator,range_rest);
            return base_generator.generate(bin_resolution, f_gen, range_first, logger); 
        }
};

template<std::size_t N, typename BaseGenerator, typename RestIntegrator>
auto regions_generator_fubini(const BaseGenerator& bg, const RestIntegrator& ri) {
    return RegionsGeneratorFubini<BaseGenerator,RestIntegrator,N>(bg,ri);
}


}