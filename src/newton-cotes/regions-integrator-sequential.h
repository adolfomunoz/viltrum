#pragma once
#include "../range.h"

namespace viltrum {

/**
 * @brief This integrator is not intended to be directly used by the user. Is a particular piece
 * common to all integrators that are based on Newton-Cotes regions. 
 */
class RegionsIntegratorSequential {

    //This is a wrapper that can be applied to used this as a regular integrator and not
    //as a region integrator. This helps reusing code, for instance on the case of
    //parallelization
    template<typename SeqRegions>
    class IntegratorWrapper {
        const RegionsIntegratorSequential& ir;
        const SeqRegions& seq_regions;
    public:
        IntegratorWrapper(const RegionsIntegratorSequential& i, const SeqRegions& seq) :
            ir(i), seq_regions(seq) {}
        
        template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
        void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
            const F& f, const Range<Float,DIM>& range, Logger& logger) const {
                return ir.integrate_regions(bins,bin_resolution,seq_regions,f,range,logger);
        }

    };

public:
    template<typename SeqRegions>
    auto integrator_wrapper(const SeqRegions& seq) const {
        return IntegratorWrapper<SeqRegions>(*this,seq);
    }

	template<typename Bins, std::size_t DIMBINS, typename SeqRegions, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate_regions(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const SeqRegions& seq_regions, const F& f, const Range<Float,DIM>& range, Logger& logger) const {

        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(bin_resolution[i]);
        std::size_t factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) factor*=bin_resolution[i];
        std::size_t progress = 0; std::size_t final_progress = seq_regions.size();
        for (const auto& r : seq_regions) {
            logger.log_progress(progress++,final_progress);
            std::array<std::size_t,DIMBINS> start_bin, end_bin;
            for (std::size_t i = 0; i<DIMBINS;++i) {
                start_bin[i] = std::size_t(std::max(Float(0),(r.range().min(i) - range.min(i))/drange[i]));
                end_bin[i]   = std::min(bin_resolution[i],std::size_t((r.range().max(i) - range.min(i))/drange[i])+1);
            }
            if (start_bin == end_bin) bins(start_bin)+=r.integral();
            else for (auto pos : multidimensional_range(start_bin, end_bin)) {
                Range<Float,DIM> binrange = range;
                for (std::size_t i=0;i<DIMBINS;++i)
                    binrange = binrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);

                Range<Float,DIM> region_bin_range = binrange.intersection(r.range());
                if (!region_bin_range.empty()) bins(pos) += double(factor)*r.integral_subrange(region_bin_range);
            } 
        }
        logger.log_progress(final_progress,final_progress);
	}
};

auto regions_integrator_sequential() {
    return RegionsIntegratorSequential();
}

}