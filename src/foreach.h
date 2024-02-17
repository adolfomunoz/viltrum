#pragma once
#include "multidimensional-range.h"
#include "log.h"
#include <thread>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <execution>

namespace viltrum {
inline struct _Sequential {} sequential;
inline struct _Parallel {} parallel;

template<std::size_t DIM, typename F>
void for_each(const _Sequential& seq, const MultidimensionalRange<DIM>& range, const F& f) {
    for (auto pos : range) f(pos);
}

template<std::size_t DIM, typename F, typename Logger>
void for_each(const _Sequential& seq, const MultidimensionalRange<DIM>& range, const F& f, Logger& logger) {
    std::size_t total = range.size();   
    std::size_t r = 0;
    for (auto pos : range) {
        logger.log_progress(r++,total);
        f(pos);
    }
    logger.log_progress(total,total);
}

template<std::size_t DIM, typename F>
void for_each(const MultidimensionalRange<DIM>& range, const F& f) {
    for_each(sequential,range,f);
}

template<std::size_t DIM, typename F, typename Logger>
void for_each(const MultidimensionalRange<DIM>& range, const F& f, Logger& logger) {
    for_each(sequential,range,f,logger);
}

template<std::size_t DIM, typename F, typename Logger>
void for_each(const _Parallel& par, const MultidimensionalRange<DIM>& range, const F& f, Logger& logger) {
    std::size_t ntasks = range.max(0)-range.min(0); 
    std::vector<std::size_t> idxs(ntasks), done(ntasks,0);
    std::iota(idxs.begin(), idxs.end(), 0);
    std::thread for_log([&done,ntasks,&logger] () {
        std::size_t i = 0;
        while (i<ntasks) {
            logger.log_progress(std::size_t(i),ntasks);
            std::this_thread::sleep_for(std::chrono::milliseconds(250));
            i = std::accumulate(done.begin(),done.end(),0);
        } 
    });
    if constexpr (DIM > 1) {
        std::array<std::size_t,DIM-1> a, b;
        for (std::size_t s = 1;s<DIM;++s) {
            a[s-1] = range.min(s);
            b[s-1] = range.max(s);
        }
        std::for_each(std::execution::par_unseq,
            idxs.begin(),idxs.end(),
            [&] (std::size_t d) {
                for (auto pos : multidimensional_range(a,b)) f(d|pos);
                ++done[d];  
            });
    } else {
        std::for_each(std::execution::par_unseq,
            idxs.begin(),idxs.end(),
            [&] (std::size_t d) { 
                f(std::array<std::size_t,1>{d}); 
                ++done[d];  
            });
    }
    for_log.join();
    logger.log_progress(ntasks,ntasks);
}

template<std::size_t DIM, typename F>
void for_each(const _Parallel& par, const MultidimensionalRange<DIM>& range, const F& f) {
    std::size_t ntasks = range.max(0)-range.min(0); 
    std::vector<std::size_t> idxs(ntasks);
    std::iota(idxs.begin(), idxs.end(), 0);
    if constexpr (DIM > 1) {
        std::array<std::size_t,DIM-1> a, b;
        for (std::size_t s = 1;s<DIM;++s) {
            a[s-1] = range.min(s);
            b[s-1] = range.max(s);
        }
        std::for_each(std::execution::par_unseq,
            idxs.begin(),idxs.end(),
            [&] (std::size_t d) {
                for (auto pos : multidimensional_range(a,b)) f(d|pos);
            });
    } else {
        std::for_each(std::execution::par_unseq,
            idxs.begin(),idxs.end(),
            [&] (std::size_t d) { 
                f(std::array<std::size_t,1>{d}); 
            });
    } 
}

}