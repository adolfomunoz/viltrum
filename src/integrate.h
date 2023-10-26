#pragma once

#include "range.h"
#include "log.h"
#include <vector>
#include <array>

namespace viltrum {

template<typename Integrator, typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
void integrate(const Integrator& integrator, Bins& bins, const std::array<std::size_t,DIMBINS>& resolution, const F& function, const Range<Float,DIM>& range, Logger& logger) {
    integrator.integrate(bins,resolution, function, range, logger);
}

template<typename Integrator, typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM>
void integrate(const Integrator& integrator, Bins& bins, const std::array<std::size_t,DIMBINS>& resolution, const F& function, const Range<Float,DIM>& range) {
    LoggerNull log;
    integrate(integrator, bins, resolution, function, range, log);
}

template<typename Integrator, typename F, typename R, typename Logger>
auto integrate(const Integrator& integrator, const F& function, const R& range, Logger& logger) {
    using T = decltype(function(range.min()));
    T sol(0);
    auto bins = [&sol] (const std::array<std::size_t,1>& i) -> T& { return sol; };
    std::array<std::size_t,1> res{1};
    integrate(integrator,bins,res,function,range,logger);
    return sol;
}

template<typename Integrator, typename F, typename R>
auto integrate(const Integrator& integrator, const F& function, const R& range) {
    LoggerNull log;
    return integrate(integrator,function,range,log);
}

template<typename Integrator, typename T, typename F, typename R, typename Logger>
void integrate(const Integrator& integrator, std::vector<T>& bins, const F& function, const R& range, Logger& logger) {
    auto b = [&bins] (const std::array<std::size_t,1>& i) -> T& { return bins[i[0]]; };
    std::array<std::size_t,1> res{bins.size()};
    integrate(integrator,b,res,function,range,logger);
}

template<typename Integrator, typename T, typename F, typename R>
void integrate(const Integrator& integrator, std::vector<T>& bins, const F& function, const R& range) {
    LoggerNull log;
    integrate(integrator,bins,function,range,log);
}




}