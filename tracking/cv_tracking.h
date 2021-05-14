#pragma once

#include "integrate-bins-adaptive.h"
#include "monte-carlo.h"
#include "sample-vector.h"

template<typename Nested, typename Error, typename ResidualStepper, typename VectorSampler>
class CV_Tracking {
	StepperAdaptive<Nested,Error> cv_stepper;
	ResidualStepper residual_stepper;
	VectorSampler vector_sampler;
	unsigned long adaptive_iterations;
    
	template<typename R,typename ResData,typename Sampler>
    struct Data {
		std::vector<R> regions;
		ResData residual_data;
		Sampler vector_sampler;
		unsigned long cv_iterations;
        Data(std::vector<R>&& rs, ResData&& rd, Sampler&& vs) : 
			regions(std::forward<std::vector<R>>(rs)),
			residual_data(std::forward<ResData>(rd)),
			vector_sampler(std::forward<Sampler>(vs)),
			cv_iterations(0) { }
    };
public:
	template<typename F, typename Float>
    auto init(const F& f, const Range<Float,1>& range) const {
		auto regions = cv_stepper.init(f,range);
        //residual_stepper.init should not do any calculation at all (should be MC)
        // return Data(std::move(regions),
		// 			residual_stepper.init(f, range),
		// 			vector_sampler(regions));
		
		auto init = residual_stepper.init(f, range);

		using VECTOR_TYPE = typename decltype(regions)::value_type;

		return Data<VECTOR_TYPE, 
					decltype(init), 
					decltype(vector_sampler(regions))>
					(std::move(regions),
					std::move(init),
					vector_sampler(regions));
    }
	
	template<typename F, typename Float, std::size_t DIM, typename R,typename ResData,typename Sampler>
    void step(const F& f, const Range<Float,DIM>& range, Data<R,ResData,Sampler>& data) const
	{
		if (data.cv_iterations<adaptive_iterations) {
			cv_stepper.step(f,range,data.regions);
			++data.cv_iterations;
		} else {
			if (data.cv_iterations == adaptive_iterations) {
				data.vector_sampler = vector_sampler(data.regions);
				++data.cv_iterations;
			}
			auto [index, probability] = data.vector_sampler.sample();
			const R& chosen_region = data.regions[index];
			residual_stepper.step([&] (const std::array<Float,DIM>& x)
		      { return (f(x) - chosen_region.approximation_at(x))/probability; },
			  chosen_region.range(), data.residual_data);
		}
    }
	
	template<typename F, typename Float, std::size_t DIM, typename R,typename ResData,typename Sampler>
    auto integral(const F& f, const Range<Float,DIM>& range, const Data<R,ResData,Sampler>& data) const {
        return cv_stepper.integral(f,range,data.regions) +
				residual_stepper.integral(f,range,data.residual_data);
    }
	
	StepperAdaptiveControlVariates(
		Nested&& nested, Error&& error, ResidualStepper&& rs, VectorSampler&& vs, unsigned long ai) :
			cv_stepper(std::forward<Nested>(nested), std::forward<Error>(error)),
			residual_stepper(std::forward<ResidualStepper>(rs)),
			vector_sampler(std::forward<VectorSampler>(vs)),
			adaptive_iterations(ai) { }
};
