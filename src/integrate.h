#pragma once

#include "range.h"

namespace viltrum {

template<typename Integrator, typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM>
void integrate(const Integrator& integrator, Bins& bins, const std::array<std::size_t,DIMBINS>& resolution, const F& function, const Range<Float,DIM>& range) {
    integrator.integrate(bins,resolution, function, range);
}

template<typename Integrator, typename F, typename Float, std::size_t DIM>
auto integrate(const Integrator& integrator, const F& function, const Range<Float,DIM>& range) {
    using T = decltype(function(range.min()));
    T sol;
    auto bins = [&sol] (const std::array<std::size_t,1>& i) -> T& { return sol; };
    std::array<std::size_t,1> res{1};
    integrate(integrator,bins,res,function,range);
    return sol;
}

}