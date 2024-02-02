#pragma once

namespace viltrum {

template<typename RegionsGenerator, typename RegionsIntegrator>
class IntegratorRegionBased {
    RegionsGenerator regions_generator;
    RegionsIntegrator regions_integrator;
public:
    IntegratorRegionBased(const RegionsGenerator& rg, const RegionsIntegrator& ri) :
        regions_generator(rg), regions_integrator(ri) {}

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {

        auto logger_generator = logger_step(logger,"region generation");
        auto seq_regions = regions_generator.generate(bin_resolution,f,range,logger_generator);
        auto logger_integration = logger_step(logger, "region integration");
        regions_integrator.integrate_regions(bins,bin_resolution,seq_regions,f,range,logger_integration);
        }
};

template<typename RegionsGenerator, typename RegionsIntegrator>
IntegratorRegionBased<RegionsGenerator,RegionsIntegrator> 
    integrator_region_based(const RegionsGenerator& rg, const RegionsIntegrator& ri) {
        return IntegratorRegionBased<RegionsGenerator,RegionsIntegrator>(rg,ri);
    }

}