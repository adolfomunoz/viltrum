#pragma once
#include <cstddef>
#include "concat.h"

/* Fubini combines different integration techniques, one of them for the first N dimensions and the other
 * for the rest of the dimensions, even if infinite.
 */

namespace viltrum {

template<typename IntegratorFirst, typename IntegratorRest, std::size_t N>
class IntegratorFubini {
    IntegratorFirst integrator_first;
    IntegratorRest integrator_rest;
public:
    IntegratorFubini(const IntegratorFirst& i_f, const IntegratorRest& i_r) :
        integrator_first(i_f), integrator_rest(i_r) {}

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {

        static_assert(N>DIMBINS,"Fubini does not work with that many dimensions on bin resolution");

        Range<Float,N> range_first;
        for (std::size_t i = 0; i<N; ++i) {
            range_first.min(i) = range.min(i);
            range_first.max(i) = range.max(i);
        }

        Range<Float,DIM-N> range_rest;
        for (std::size_t i = N; i<DIM; ++i) {
            range_rest.min(i-N) = range.min(i);
            range_rest.max(i-N) = range.max(i);
        }      
        viltrum::integrate(integrator_first,bins,bin_resolution,
            [&] (const std::array<float,N>& x) {
                return integrate(integrator_rest,[&f,&x] (const std::array<float,DIM-N>& xr) {
                    return f(x | xr);
                }, range_rest);
            }, range_first, logger);
	}

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const RangeInfinite<Float>& range, Logger& logger) const {

        static_assert(N>DIMBINS,"Fubini does not work with that many dimensions on bin resolution");

        std::array<Float,N> first_min, first_max;
        for (std::size_t i = 0; i<N; ++i) {
            first_min[i] = range.min(i);
            first_max[i] = range.max(i);
        }

        std::vector<Float> rest_min = range.min();
        std::vector<Float> rest_max = range.min();
        for (std::size_t i = 0; i<N; ++i) {
            if (!rest_min.empty()) rest_min.erase(rest_min.begin());
            if (!rest_max.empty()) rest_max.erase(rest_max.begin());
        }
   
        viltrum::integrate(integrator_first,bins,bin_resolution,
            [&] (const std::array<Float,N>& x) {
                return viltrum::integrate(integrator_rest,[&f,&x] (const auto& xr) {
                    return f(concat(x,xr));
                }, range_infinite(rest_min,rest_max));
            }, Range<Float,N>(first_min,first_max), logger);
    }
};

template<std::size_t N, typename IntegratorFirst, typename IntegratorRest>
IntegratorFubini<IntegratorFirst,IntegratorRest,N> integrator_fubini(IntegratorFirst&& ifirst, IntegratorRest&& irest) {
    return IntegratorFubini<std::decay_t<IntegratorFirst>,std::decay_t<IntegratorRest>,N>
        (std::forward<IntegratorFirst>(ifirst),std::forward<IntegratorRest>(irest));
}

};