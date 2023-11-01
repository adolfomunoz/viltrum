#pragma once

#include "range.h"
#include "log.h"
#include <vector>
#include <array>

namespace viltrum {

namespace detail {
    template<typename P, typename F, typename Enable = void>
    struct Integrand {
        static const F& adapt(const F& f) { return f; } 
    };
  
    template<typename P, typename F>
    struct Integrand<P,F,typename  std::enable_if_t<std::is_invocable_v<F,P>> > {
        static auto adapt(const F& f) { 
            return [&f] (const std::array<P,1>& x) { return f(x[0]);};   
        } 
    };

    template<typename P, typename F>
    struct Integrand<P,F,typename  std::enable_if_t<std::is_invocable_v<F,P,P>> > {
        static auto adapt(const F& f) { 
            return [&f] (const std::array<P,2>& x) { return f(x[0],x[1]);};   
        } 
    };

    template<typename P, typename F>
    struct Integrand<P,F,typename  std::enable_if_t<std::is_invocable_v<F,P,P,P>> > {
        static auto adapt(const F& f) { 
            return [&f] (const std::array<P,3>& x) { return f(x[0],x[1],x[2]);};   
        } 
    };

    template<typename P, typename F>
    struct Integrand<P,F,typename  std::enable_if_t<std::is_invocable_v<F,P,P,P,P>> > {
        static auto adapt(const F& f) { 
            return [&f] (const std::array<P,4>& x) { return f(x[0],x[1],x[2],x[3]);};   
        } 
    };

    template<typename P, typename F>
    auto adapt(const F& f) {
        return Integrand<P,std::decay_t<F>>::adapt(f);
    }   
} 

template<typename Integrator, typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
void integrate(const Integrator& integrator, Bins& bins, const std::array<std::size_t,DIMBINS>& resolution, 
            const F& function, const Range<Float,DIM>& range, Logger& logger,
            std::enable_if_t<std::is_convertible_v<
                std::decay_t<decltype(detail::adapt<Float>(function)(std::array<Float,DIM>()))>,
                std::decay_t<decltype(bins(std::array<std::size_t,DIMBINS>()))>>,int> dummy=0) {
    integrator.integrate(bins,resolution, detail::adapt<Float>(function), range, logger);
}



template<typename Integrator, typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM>
void integrate(const Integrator& integrator, Bins& bins, const std::array<std::size_t,DIMBINS>& resolution, const F& function, const Range<Float,DIM>& range) {
    LoggerNull log;
    integrate(integrator, bins, resolution, function, range, log);
}

template<typename Integrator, typename F, typename R, typename Logger>
auto integrate(const Integrator& integrator, const F& function, const R& range, Logger& logger) {
    using T = decltype(detail::adapt<typename R::value_type>(function)(range.min()));
    T sol(0.0);
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