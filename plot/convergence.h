#pragma once

#include <chrono>
#include <functional>
#include "../quadrature/integrate-bins.h"

namespace {
template<typename F,std::size_t DIM>
class FunctionWrapper {
	F f;
	mutable std::shared_ptr<unsigned long> nsamples;
public:
	FunctionWrapper(const F& f) : 
		f(f), nsamples(std::make_shared<unsigned long>(0)) {}
	FunctionWrapper(F&& f) : 
		f(std::forward<F>(f)), nsamples(std::make_shared<unsigned long>(0)) {}
	
	double operator()(const std::array<double,DIM>& pos) const {
		double value = std::apply(f,pos);
		++(*nsamples);
		return value;
	}
	
	float operator()(const std::array<float,DIM>& pos) const {
		float value = std::apply(f,pos);
		++(*nsamples);
		return value;
	}	
	unsigned long samples() const { return *nsamples; }
};
}

template<typename Stepper, typename F, typename Error, std::size_t DIM>
svg_cpp_plot::_2d::polyline convergence_samples(const Stepper& stepper, const F& func, const Range<double,DIM>& range, const Error& error, unsigned long max_samples, std::size_t averaging = 1, unsigned int resolution = 100) {
	FunctionWrapper<F,DIM> f(func); unsigned long samples_plot = 0;
	svg_cpp_plot::_2d::polyline samples_vs_error;
    using StepperData = decltype(stepper.init(f,range));
    std::vector<StepperData> data;
    for (std::size_t i = 0; i<averaging; ++i) data.push_back(stepper.init(f,range));
    double err;
    while (samples_plot <= max_samples) {
        for (auto& d : data) stepper.step(f,range,d);
        if (f.samples() > samples_plot*averaging) {
            err = 0; for (auto d : data) err += (error(stepper.integral(f,range,d))/double(averaging)); 
            samples_vs_error.add_point(std::log(f.samples()/averaging),std::log(err));
            samples_plot += max_samples/resolution;
        }
    }
    return samples_vs_error;
}



template<typename Stepper, typename F, typename Error, std::size_t DIM>
svg_cpp_plot::_2d::polyline convergence_time(const Stepper& stepper, const F& func, const Range<double,DIM>& range, const Error& error, double max_time, std::size_t averaging = 1, unsigned int resolution = 100) {
    double err;
	svg_cpp_plot::_2d::polyline time_vs_error;
    using StepperData = decltype(stepper.init(FunctionWrapper<F,DIM>(func),range));
    for (double time_plot = 0; time_plot < max_time; time_plot += max_time/resolution) {
	    FunctionWrapper<F,DIM> f(func);
	    auto start = std::chrono::steady_clock::now();
        std::vector<StepperData> data;
        for (std::size_t i = 0; i<averaging; ++i) data.push_back(stepper.init(f,range));
        for (auto& d : data) stepper.step(f,range,d);
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
        while (elapsed.count() <= time_plot*averaging) {
            for (auto& d : data) stepper.step(f,range,d);
            end = std::chrono::steady_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
        }
        err = 0; for (auto d : data) err += (error(stepper.integral(f,range,d))/double(averaging)); 
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
        time_vs_error.add_point(std::log(elapsed.count()/averaging),std::log(err));
    }
    return time_vs_error;
}

template<typename Stepper, typename F, typename Error, std::size_t DIM>
svg_cpp_plot::_2d::polyline convergence(bool samples, const Stepper& stepper, const F& func, const Range<double,DIM>& range, const Error& error, double max_what, std::size_t averaging = 1, unsigned int resolution = 100) {
	if (samples) 
		return convergence_samples(stepper,func,range,error,(unsigned long)(max_what),averaging,resolution);
	else
		return convergence_time(stepper,func,range,error,max_what,averaging,resolution);
		
}


template<typename Stepper, typename F, typename Error, std::size_t DIM>
svg_cpp_plot::_2d::polyline convergence_bins_samples(const Stepper& stepper, const F& func, const Range<double,DIM>& range, const Error& error, unsigned long max_samples, unsigned long bins, std::size_t averaging = 1, unsigned int resolution = 100) {
    std::vector<double> sol(bins); unsigned long samples_plot = 0;
	FunctionWrapper<F,DIM> f(func);
	svg_cpp_plot::_2d::polyline samples_vs_error;
    using StepperData = decltype(stepper.init(vector_resolution(sol),f,range));
    std::vector<StepperData> data;
    for (std::size_t i = 0; i<averaging; ++i) data.push_back(stepper.init(vector_resolution(sol),f,range));
    double err;
    while (samples_plot <= max_samples) {
        for (auto& d : data) stepper.step(vector_resolution(sol),f,range,d);
        if (f.samples() > samples_plot*averaging) {
            err = 0; for (auto d : data) {
                std::fill(sol.begin(),sol.end(),0.0);
                auto bins = vector_bins(sol);
                stepper.integral(bins, vector_resolution(sol), f, range, d);
                err += (error(sol)/double(averaging));
            } 
            samples_vs_error.add_point(std::log(std::max(1.e-20,double(f.samples())/averaging)),std::log(std::max(1.e-20,err)));
            samples_plot += max_samples/resolution;
        }
    }
    return samples_vs_error;
}

template<typename IntegratorSpp, typename F, typename Error, std::size_t DIM>
svg_cpp_plot::_2d::polyline convergence_bins_integrator_samples(const IntegratorSpp& integrator_spp, const F& func, const Range<double,DIM>& range, const Error& error, unsigned long max_samples, unsigned long bins, std::size_t averaging = 1, unsigned int resolution = 100) {
    std::vector<double> sol(bins); double samples_plot;
    double iteration_factor = std::pow(double(max_samples),1.0/double(resolution));
    auto vbins = vector_bins(sol);
	FunctionWrapper<F,DIM> f(func); 
    svg_cpp_plot::_2d::polyline samples_vs_error;
    double err;
    for (samples_plot = 1; samples_plot <= max_samples; samples_plot *= iteration_factor) {
        err = 0;
        auto integrator = integrator_spp((unsigned long)(samples_plot));
        for (std::size_t a = 0; a<averaging; ++a) {
            integrator.integrate(vbins,vector_resolution(sol), f, range);
            err += error(sol);
        }
        samples_vs_error.add_point(std::log(samples_plot),std::log(err/averaging));
    }
    return samples_vs_error;
}


template<typename Stepper, typename F, typename Error, std::size_t DIM>
svg_cpp_plot::_2d::polyline convergence_bins_time(const Stepper& stepper, const F& func, const Range<double,DIM>& range, const Error& error, double max_time, unsigned long bins, std::size_t averaging = 1, unsigned int resolution = 100) {
    std::vector<double> sol(bins);
    double err;
	svg_cpp_plot::_2d::polyline time_vs_error;
    using StepperData = decltype(stepper.init(vector_resolution(sol),FunctionWrapper<F,DIM>(func),range));
    for (double time_plot = 0; time_plot < max_time; time_plot += max_time/resolution) {
	    FunctionWrapper<F,DIM> f(func);
	    auto start = std::chrono::steady_clock::now();
        std::vector<StepperData> data;
        for (std::size_t i = 0; i<averaging; ++i) data.push_back(stepper.init(vector_resolution(sol),f,range));
        for (auto& d : data) stepper.step(vector_resolution(sol),f,range,d);
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
        while (elapsed.count() <= time_plot*averaging) {
            for (auto& d : data) stepper.step(vector_resolution(sol),f,range,d);
            end = std::chrono::steady_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
        }
        err = 0; for (auto d : data) {
            std::fill(sol.begin(),sol.end(),0.0);
            auto bins = vector_bins(sol);
            stepper.integral(bins, vector_resolution(sol), f, range, d);
            err += (error(sol)/double(averaging));
        } 
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
        time_vs_error.add_point(std::log(elapsed.count()/averaging),std::log(err));
    }
    return time_vs_error;
}

template<typename Stepper, typename F, typename Error, std::size_t DIM>
svg_cpp_plot::_2d::polyline convergence_bins(bool samples, const Stepper& stepper, const F& func, const Range<double,DIM>& range, const Error& error, double max_what, unsigned long bins, std::size_t averaging = 1, unsigned int resolution = 100) {
	if (samples) 
		return convergence_bins_samples(stepper,func,range,error,(unsigned long)(max_what),bins,averaging,resolution);
	else
		return convergence_bins_time(stepper,func,range,error,max_what,bins,averaging,resolution);
	
}


