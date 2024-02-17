#pragma once
#include <cstddef>
#include <array>
#include <vector>
#include "sequence_from_iterators.h"


/* Reorder reorders the dimensions of integration, so it is well combined with Fubini when needing
 * to integrate some particular dimensions with one technique, without being the first ones
 */

namespace viltrum {

template<typename Integrator, std::size_t N>
class IntegratorReorder {
    Integrator integrator;
    std::array<std::size_t,N> first_dimensions;
public:
    IntegratorReorder(const Integrator& i, const std::array<std::size_t,N>& fd) :
        integrator(i), first_dimensions(fd) {}

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {

        std::array<Float,DIM> reorder_min = range.min(); 
        std::array<Float,DIM> reorder_max = range.max();
        for (std::size_t i = 0; i<N; ++i) {
            if (first_dimensions[i]<DIM) {
                std::swap(reorder_min[i],reorder_min[first_dimensions[i]]);
                std::swap(reorder_max[i],reorder_max[first_dimensions[i]]);
            }
        }
        Range<Float,DIM> range_reordered(reorder_min,reorder_max);


        const auto f_reordered = [&] (const std::array<Float,DIM>& x) {
            std::array<Float,DIM> x_reordered = x;
            for (std::size_t i = 0; i<N; ++i) {
                if (first_dimensions[i]<DIM) {
                    std::swap(x_reordered[i],x_reordered[first_dimensions[i]]);
                }
            }
            return f(x_reordered);
        };

        integrator.integrate(bins,bin_resolution, f_reordered, range_reordered, logger);
	}

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const RangeInfinite<Float>& range, Logger& logger) const {

        std::size_t size_reorder = 0; //For checking max size of the reorder vector
        for (std::size_t i : first_dimensions) if (size_reorder<=i) size_reorder = (i+1);
        
        std::vector<Float> reorder_min = range.min(); 
        if (reorder_min.size()<size_reorder) reorder_min.resize(size_reorder,0);
        std::vector<Float> reorder_max = range.max();
        if (reorder_max.size()<size_reorder) reorder_max.resize(size_reorder,1);

        for (std::size_t i = 0; i<N; ++i) {
            std::swap(reorder_min[i],reorder_min[first_dimensions[i]]);
            std::swap(reorder_max[i],reorder_max[first_dimensions[i]]);
        } 
        RangeInfinite<Float> range_reordered(reorder_min,reorder_max);       
        
        const auto f_reordered = [&] (const auto& x) {
            std::vector<Float> x_reordered(size_reorder);
            auto it = x.begin();
            for (Float& xr : x_reordered) { xr = (*it); ++it; }
            
            for (std::size_t i = 0; i<N; ++i) {
                std::swap(x_reordered[i],x_reordered[first_dimensions[i]]);
            }
            return f(concat(x_reordered,sequence_from_iterators(it,x.end())));
        };

        integrator.integrate(bins,bin_resolution, f_reordered, range_reordered, logger);
    }
};

template<typename Integrator, std::size_t N>
IntegratorReorder<Integrator,N> integrator_reorder(Integrator&& i, const std::array<std::size_t,N>& fd) {
    return IntegratorReorder<std::decay_t<Integrator>,N>(std::forward<Integrator>(i),fd);
}

};