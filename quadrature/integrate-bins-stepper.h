#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <chrono>
#include "integrate-bins.h"
#include "vector-dimensions.h"

namespace viltrum {

template<typename Stepper>
class IntegratorBinsStepper {
    Stepper stepper;
    unsigned long iterations;
public:

	template<typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM>
	void integrate(Bins& bins, const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range) const {
        auto data = stepper.init(bin_resolution, f, range);
        for (unsigned long i = 0; i<iterations;++i) stepper.step(bin_resolution,f,range,data);
        stepper.integral(bins,bin_resolution,f,range,data);
    }

    IntegratorBinsStepper(Stepper&& s, unsigned long i) :
        stepper(std::forward<Stepper>(s)), iterations(i) { }
};

template<typename Stepper>
auto integrator_bins_stepper(Stepper&& stepper, unsigned long iterations) {
    return IntegratorBinsStepper<std::decay_t<Stepper>>(std::forward<Stepper>(stepper), iterations);
}

template<typename StepperPerBin>
class StepperBinsPerBin {
    StepperPerBin bin_stepper;
public:
    template<std::size_t DIMBINS, typename F, typename Float, std::size_t DIM>
    auto init(const std::array<std::size_t,DIMBINS>& resolution, const F& f, const Range<Float,DIM>& range) const {
        using StepperData = decltype(bin_stepper.init(f,range));
        vector_dimensions<StepperData,DIMBINS> data(resolution);
        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(resolution[i]);
        
        for (auto pos : multidimensional_range(resolution)) {
            Range<Float,DIM> subrange = range;
            for (std::size_t i=0;i<DIMBINS;++i)
                subrange = subrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);
            data[pos] = bin_stepper.init(f,subrange);
        }
        return data;
    }

    template<std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename StepperData>
    void step(const std::array<std::size_t,DIMBINS>& resolution, const F& f, const Range<Float,DIM>& range, vector_dimensions<StepperData,DIMBINS>& data) const {
        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(resolution[i]);
        
        for (auto pos : multidimensional_range(resolution)) {
            Range<Float,DIM> subrange = range;
            for (std::size_t i=0;i<DIMBINS;++i) 
                subrange = subrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);
            bin_stepper.step(f,subrange,data(pos));
        }
    }

    template<typename Bins, typename F, typename Float, std::size_t DIM, std::size_t DIMBINS, typename StepperData>
    void integral(Bins& bins, const std::array<std::size_t,DIMBINS>& resolution, const F& f, const Range<Float,DIM>& range, const vector_dimensions<StepperData,DIMBINS>& data) const {
        std::array<Float,DIMBINS> drange;
        for (std::size_t i=0;i<DIMBINS;++i) drange[i] = (range.max(i) - range.min(i))/Float(resolution[i]);
        for (auto pos : multidimensional_range(resolution)) {
            Range<Float,DIM> subrange = range;
            for (std::size_t i=0;i<DIMBINS;++i) 
                subrange = subrange.subrange_dimension(i,range.min(i)+pos[i]*drange[i],range.min(i)+(pos[i]+1)*drange[i]);
            bins(pos) = data.size()*bin_stepper.integral(f,subrange,data(pos));
        }
    }

    StepperBinsPerBin(StepperPerBin&& bs) : bin_stepper(std::forward<StepperPerBin>(bs)) { }
};

template<typename StepperPerBin>
auto stepper_bins_per_bin(StepperPerBin&& bin_stepper) {
    return StepperBinsPerBin<std::decay_t<StepperPerBin>>(std::forward<StepperPerBin>(bin_stepper));
}

template<typename Stepper, typename Bins, std::size_t DIMBINS, typename F, typename Float, std::size_t DIM>
void integrate_bins_stepper_progression(const std::string& name, const Stepper& stepper, unsigned long iterations, Bins& bins, const std::array<std::size_t,DIMBINS>& resolution, const F& f, const Range<Float,DIM>& range) {
    std::cerr<<name<<" - \r";
	auto start = std::chrono::steady_clock::now();
    auto data = stepper.init(resolution, f, range);
    for (unsigned long i = 0; i<iterations;++i) {
        stepper.step(resolution,f,range,data);
	    auto end = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(end - start); 
        std::cerr<<name<<" - \t"<<std::fixed<<std::setprecision(2)<<std::setw(6)<<(100.0f*float(i)/float(iterations))<<"%\t("<<std::setprecision(3)<<std::setw(6)<<elapsed.count()<<" seconds)\r";
    }   
    stepper.integral(bins,resolution,f,range,data);
    auto end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(end - start); 
    std::cout<<name<<" - \t[DONE] \t("<<std::setprecision(3)<<std::setw(6)<<elapsed.count()<<" seconds)"<<std::endl;
}

}


