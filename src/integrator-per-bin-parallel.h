#pragma once
#include "range.h"
#include "multidimensional-range.h"
#include "integrate.h"
#include <thread>
#include <atomic>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <execution>

namespace viltrum {

template<typename Integrator>
class IntegratorPerBinParallel {
    Integrator bin_integrator;

public:
	template<typename Bins, std::size_t DIMBINS, typename F, typename R, typename Logger>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const R& range, Logger& logger) const {
        std::size_t factor(1);   
        for (std::size_t i = 0; i < DIMBINS; ++i) factor*=bin_resolution[i];
        std::array<typename R::value_type,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/((typename R::value_type)(bin_resolution[i]));

        std::size_t ntasks = bin_resolution[0]; 
        std::vector<std::size_t> idxs(bin_resolution[0]), done(bin_resolution[0],0);
        std::iota(idxs.begin(), idxs.end(), 0);
//        std::atomic<std::size_t> i = 0;
        std::thread for_log([&done,ntasks,&logger] () {
            std::size_t i = 0;
            while (i<ntasks) {
                logger.log_progress(std::size_t(i),ntasks);
                std::this_thread::sleep_for(std::chrono::milliseconds(250));
                i = std::accumulate(done.begin(),done.end(),0);
            } 
        });
        if constexpr (DIMBINS > 1) {
            std::array<std::size_t,DIMBINS-1> task_dims;
            for (std::size_t s = 1;s<DIMBINS;++s)
                task_dims[s-1] = bin_resolution[s];
            std::for_each(std::execution::par_unseq,
                idxs.begin(),idxs.end(),
                [&] (std::size_t d) {
                    for (auto pos : multidimensional_range(task_dims)) {
                        R subrange = range;
                        subrange=range.subrange_dimension(0,range.min(0)+d*drange[0],range.min(0)+(d+1)*drange[0]);
                        for (std::size_t n=1;n<DIMBINS;++n)
                            subrange = subrange.subrange_dimension(n,range.min(n)+pos[n-1]*drange[n],range.min(n)+(pos[n-1]+1)*drange[n]);
                        bins(d|pos) = double(factor)*viltrum::integrate(bin_integrator,f,subrange);
                    }
                    ++done[d];         
                });  
        } else {
            std::for_each(std::execution::par_unseq,
                idxs.begin(),idxs.end(),
                [&] (std::size_t d) {
                    R subrange = range;
                    subrange=range.subrange_dimension(0,range.min(0)+d*drange[0],range.min(0)+(d+1)*drange[0]);
                    bins(std::array<std::size_t,1>{d}) = double(factor)*viltrum::integrate(bin_integrator,f,subrange);
                    ++done[d];         
                });  
        } 
        for_log.join();  
        logger.log_progress(factor,factor);
	}

    IntegratorPerBinParallel(Integrator&& pi) : 
	    bin_integrator(std::forward<Integrator>(pi)) { }
    IntegratorPerBinParallel(const Integrator& pi) : 
	    bin_integrator(pi) { }

};

template<typename Integrator>
auto integrator_per_bin_parallel(Integrator&& i) {
    return IntegratorPerBinParallel<std::decay_t<Integrator>>(std::forward<Integrator>(i));
}

}