#pragma once
#include <cstddef>
#include "concat.h"
#include "reorder.h"
#include <tuple>
#include "../range.h"
#include "../range-infinite.h"
#include "../array.h"



/* Fubini combines different integration techniques, one of them for the first N dimensions and the other
 * for the rest of the dimensions, even if infinite.
 */

namespace viltrum {

template<std::size_t N, typename Float, std::size_t DIM>
std::tuple<Range<Float,N>, Range<Float,DIM-N>> range_split_at(const Range<Float,DIM>& range) {
    std::array<Float,N> range_first_min, range_first_max;
    for (std::size_t i = 0; i<N; ++i) {
        range_first_min[i] = range.min(i);
        range_first_max[i] = range.max(i);
    }
    std::array<Float,DIM-N> range_rest_min, range_rest_max;
    for (std::size_t i = N; i<DIM; ++i) {
        range_rest_min[i-N] = range.min(i);
        range_rest_max[i-N] = range.max(i);
    }      
    return std::tuple(Range<Float,N>(range_first_min, range_first_max),
                      Range<Float,DIM-N>(range_rest_min, range_rest_max));
}

template<std::size_t N, typename Float>
std::tuple<Range<Float,N>, RangeInfinite<Float>> range_split_at(const RangeInfinite<Float>& range) {
    std::array<Float,N> range_first_min, range_first_max;
    for (std::size_t i = 0; i<N; ++i) {
        range_first_min[i] = range.min(i);
        range_first_max[i] = range.max(i);
    }
        std::vector<Float> range_rest_min = range.min();
        std::vector<Float> range_rest_max = range.max();
        for (std::size_t i = 0; i<N; ++i) {
            if (!range_rest_min.empty()) range_rest_min.erase(range_rest_min.begin());
            if (!range_rest_max.empty()) range_rest_max.erase(range_rest_max.begin());
        }    
        return std::tuple(Range<Float,N>(range_first_min, range_first_max),
                     RangeInfinite<Float>(range_rest_min, range_rest_max));
}

template<std::size_t N, typename F, typename Integrator, typename Float>
auto function_split_and_integrate_at(const F& f, const Integrator& integrator, const Range<Float,0>& range) {
    return f; 
} 

template<std::size_t N, typename F, typename Integrator, typename Float, std::size_t DIM>
auto function_split_and_integrate_at(const F& f, const Integrator& integrator, const Range<Float,DIM>& range) {
    return [=] (const std::array<Float,N>& x) {
        return viltrum::integrate(integrator,[&f,&x] (const std::array<float,DIM>& xr) {
            return f( (x | xr) );
        }, range);
    };  
}

template<std::size_t N, typename F, typename Integrator, typename Float>
auto function_split_and_integrate_at(const F& f, const Integrator& integrator, const RangeInfinite<Float>& range) {
    return [=] (const std::array<Float,N>& x) {
        return viltrum::integrate(integrator,[&f,&x] (const auto& xr) {
            return f(concat(x,xr));
        }, range);
    }; 
} 


template<typename IntegratorFirst, typename IntegratorRest, std::size_t N>
class IntegratorFubini {
    IntegratorFirst integrator_first;
    IntegratorRest integrator_rest;
public:
    IntegratorFubini(const IntegratorFirst& i_f, const IntegratorRest& i_r) :
        integrator_first(i_f), integrator_rest(i_r) {}

	template<typename Bins, std::size_t DIMBINS, typename F, typename Range, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range& range, Logger& logger) const {

        static_assert(N>=DIMBINS,"Fubini does not work with that many dimensions on bin resolution");

        auto [range_first, range_rest] = range_split_at<N>(range);

        viltrum::integrate(integrator_first,bins,bin_resolution,
            function_split_and_integrate_at<N>(f,integrator_rest,range_rest),
            range_first, logger);
	}
};

template<std::size_t N, typename IntegratorFirst, typename IntegratorRest>
IntegratorFubini<IntegratorFirst,IntegratorRest,N> integrator_fubini(IntegratorFirst&& ifirst, IntegratorRest&& irest) {
    return IntegratorFubini<std::decay_t<IntegratorFirst>,std::decay_t<IntegratorRest>,N>
        (std::forward<IntegratorFirst>(ifirst),std::forward<IntegratorRest>(irest));
}

template<typename IntegratorFirst, typename IntegratorRest>
auto integrator_fubini(const std::tuple<std::size_t>& dims, IntegratorFirst&& ifirst, IntegratorRest&& irest) {
    std::vector<std::tuple<size_t,size_t>> reorder{ {0, std::get<0>(dims)}};
    return integrator_reorder(integrator_fubini<1>(std::forward<IntegratorFirst>(ifirst), std::forward<IntegratorRest>(irest)),{0,std::get<0>(dims)});
} 

template<typename IntegratorFirst, typename IntegratorRest>
auto integrator_fubini(const std::tuple<std::size_t,std::size_t>& dims, IntegratorFirst&& ifirst, IntegratorRest&& irest) {
    std::vector<std::tuple<size_t,size_t>> reorder{ {0, std::get<0>(dims)}, {1, std::get<1>(dims)}};
    return integrator_reorder(integrator_fubini<2>(std::forward<IntegratorFirst>(ifirst), std::forward<IntegratorRest>(irest)),reorder);
} 

template<typename IntegratorFirst, typename IntegratorRest>
auto integrator_fubini(const std::tuple<std::size_t,std::size_t,std::size_t>& dims, IntegratorFirst&& ifirst, IntegratorRest&& irest) {
    std::vector<std::tuple<size_t,size_t>> reorder{ {0, std::get<0>(dims)}, {1, std::get<1>(dims)}, {2, std::get<2>(dims)}};
    return integrator_reorder(integrator_fubini<3>(std::forward<IntegratorFirst>(ifirst), std::forward<IntegratorRest>(irest)),reorder);
} 

template<typename IntegratorFirst, typename IntegratorRest>
auto integrator_fubini(const std::tuple<std::size_t,std::size_t,std::size_t,std::size_t>& dims, IntegratorFirst&& ifirst, IntegratorRest&& irest) {
    std::vector<std::tuple<size_t,size_t>> reorder{ {0, std::get<0>(dims)}, {1, std::get<1>(dims)}, {2, std::get<2>(dims)}, {3, std::get<3>(dims)}};
    return integrator_reorder(integrator_fubini<4>(std::forward<IntegratorFirst>(ifirst), std::forward<IntegratorRest>(irest)),reorder);
} 

};