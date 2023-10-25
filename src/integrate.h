#pragma once

#include "range.h"
#include <vector>
#include <array>

namespace viltrum {

template<typename Integrator, typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM>
void integrate(const Integrator& integrator, Bins& bins, const std::array<std::size_t,DIMBINS>& resolution, const F& function, const Range<Float,DIM>& range) {
    integrator.integrate(bins,resolution, function, range);
}

template<typename Integrator, typename F, typename R>
auto integrate(const Integrator& integrator, const F& function, const R& range) {
    using T = decltype(function(range.min()));
    T sol;
    auto bins = [&sol] (const std::array<std::size_t,1>& i) -> T& { return sol; };
    std::array<std::size_t,1> res{1};
    integrate(integrator,bins,res,function,range);
    return sol;
}

template<typename Integrator, typename T, typename F, typename R>
void integrate(const Integrator& integrator, std::vector<T>& bins, const F& function, const R& range) {
    auto b = [&bins] (const std::array<std::size_t,1>& i) -> T& { return bins[i[0]]; };
    std::array<std::size_t,1> res{bins.size()};
    integrate(integrator,b,res,function,range);
}



}