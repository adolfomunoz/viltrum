#pragma once
#include "../range.h"
#include "../integrator-per-bin-parallel.h"

namespace viltrum {

/**
 * @brief This integrator is not intended to be directly used by the user. Is a particular piece
 * common to all integrators that are based on Newton-Cotes regions. The integrand is therefore ignored
 * 
 * @tparam SeqRegions Any STL-like linear sequence of regions
 */
template<typename SeqRegions>
class IntegratorRegions {
    SeqRegions seq_regions;
public:
    IntegratorRegions(const SeqRegions& sr) : seq_regions(sr) {}
    IntegratorRegions(SeqRegions&& sr) : seq_regions(std::forward<SeqRegions>(sr)) {}

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {
        //Note that this function is already ignored, we expect the input to be the regions which
        // have already been precalculated
        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(bin_resolution[i]);
        std::size_t factor = 1;
        for (std::size_t i=0;i<DIMBINS;++i) factor*=bin_resolution[i];
        std::size_t progress = 0; std::size_t final_progress = seq_regions.size();

        for (const auto& r : seq_regions) {
            logger.log_progress(progress++,final_progress);
            std::array<std::size_t,DIMBINS> start_bin, end_bin;
            for (std::size_t i = 0; i<DIMBINS;++i) {
                start_bin[i] = std::max(std::size_t(0),std::size_t((r.range().min(i) - range.min(i))/drange[i]));
                end_bin[i]   = std::min(bin_resolution[i],std::size_t((r.range().max(i) - range.min(i))/drange[i])+1);
            }
            if (start_bin == end_bin) bins(start_bin)+=r.integral();
            else for (auto pos : multidimensional_range(start_bin, end_bin)) {
                auto subrange = range;
                for (std::size_t i=0;i<DIMBINS;++i)
                    subrange = subrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);
                bins(pos) += double(factor)*r.integral_subrange(subrange);
            } 
        }
        logger.log_progress(final_progress,final_progress);
	}
};

template<typename RegionSeq, typename Bins, std::size_t DIMBINS, typename Float, std::size_t DIM, typename Logger>
void integrate_regions(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const RegionSeq& region_seq, const Range<Float,DIM>& range, bool parallel, Logger& logger) {

    auto dummy_function = 
        [] (const std::array<Float,DIM>& x) { return typename std::decay_t<decltype(region_seq.front())>::value_type(); };
    if (parallel) {
        integrate(integrator_per_bin_parallel(IntegratorRegions<RegionSeq>(region_seq)),
            bins,bin_resolution,dummy_function,range,logger);
    } else {
        integrate(IntegratorRegions<RegionSeq>(region_seq),
            bins,bin_resolution,dummy_function,range,logger);
    }
}

}