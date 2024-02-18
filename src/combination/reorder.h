#pragma once
#include <cstddef>
#include <array>
#include <vector>
#include "sequence_from_iterators.h"


/* Reorder reorders the dimensions of integration, so it is well combined with Fubini when needing
 * to integrate some particular dimensions with one technique, without being the first ones
 */

namespace viltrum {

template<typename Integrator>
class IntegratorReorder {
    Integrator integrator;
    std::vector<std::tuple<std::size_t,std::size_t>> to_from;

    template<typename V>
    std::vector<typename V::value_type> reorder(const V& v, const typename V::value_type& def = typename V::value_type()) const {

        std::vector<bool> tos(v.size(),false), froms(v.size(),false);
        for (auto [to, from] : to_from) {
            if (from<v.size()) { 
                if (to>=tos.size()) tos.resize(to+1,false);
                froms[from]=true; tos[to]=true;
            } 
        }
        std::size_t ntos = 0, nfroms = 0;
        std::for_each(tos.begin(),tos.end(),[&ntos] (bool b) { if (b) ++ntos; });
        std::for_each(froms.begin(),froms.end(),[&nfroms] (bool b) { if (b) ++nfroms; });

        std::vector<typename V::value_type> sol(std::max(tos.size(),ntos + v.size()-nfroms),def);
        for (auto [to, from] : to_from) {
            if ((from<v.size()) && (to<sol.size())) {
                sol[to]=v[from];    
            } 
        }
        std::size_t from, to = 0;
        for (from = 0; from<v.size(); ++from) {
            if (!froms[from]) { //We don't have this value yet
                while ((to<sol.size()) && tos[to]) ++to; //Find where to put it
                if (to<sol.size()) sol[to++] = v[from]; //Actually put it  
            }   
        } 
        return sol;   
    }

public:
    IntegratorReorder(const Integrator& i, const std::vector<std::tuple<std::size_t,std::size_t>>& fd) :
        integrator(i), to_from(fd) {}

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {

        auto range_reordered = viltrum::range<DIM>(reorder(range.min(),Float(0)),reorder(range.max(),Float(1)));
        const auto f_reordered = [&] (const std::array<Float,DIM>& x) {
            auto re = reorder(x);
            std::array<Float,DIM> x_re;
            for (std::size_t i = 0; (i<DIM) && (i<re.size()); ++i) x_re[i] = re[i];  
            return (f(x_re));
        };

        integrator.integrate(bins,bin_resolution, f_reordered, range_reordered, logger);
	}

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const RangeInfinite<Float>& range, Logger& logger) const {

        RangeInfinite<Float> range_reordered(reorder(range.min(),Float(0)),reorder(range.max(),Float(1)));       
         std::size_t size_reorder = 0;
        for (auto [to, from] : to_from) {
            if (size_reorder<=from) {
                size_reorder = from + 1; 
            } 
        }       
        const auto f_reordered = [&] (const auto& x) {
            std::vector<Float> xr(size_reorder);
            auto it = x.begin();
            for (Float& v : xr) { v = (*it); ++it; }
/*           for (Float v : xr) std::cerr<<v<<" ";
            std::cerr<<" -> ";
            for (Float v : reorder(xr)) std::cerr<<v<<" ";
            std::cerr<<std::endl; */
            return f(concat(reorder(xr),sequence_from_iterators(it,x.end())));
        };

        integrator.integrate(bins,bin_resolution, f_reordered, range_reordered, logger);
    }
};

template<typename Integrator>
IntegratorReorder<Integrator> integrator_reorder(Integrator&& i, const std::vector<std::tuple<std::size_t,std::size_t>>& fd) {
    return IntegratorReorder<std::decay_t<Integrator>>(std::forward<Integrator>(i),fd);
}

};