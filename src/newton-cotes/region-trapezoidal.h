#pragma once

#include <array>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include "../multiarray/multiarray.h"
#include "../multiarray/fill.h"
#include "../range.h"

#ifndef VILTRUM_MAX_DIMENSIONS_REGION
#define VILTRUM_MAX_DIMENSIONS_REGION 6
#endif

namespace viltrum {

/**
 * Region that is always approximated with the trapezoidal rule: the function is the multilinear
 * interpolation of its values at the 2^DIM corners of the region.
 *
 * Every operation is a weighted sum of the corner values, where the weight of a corner is the 
 * product of one-dimensional weights (one pair [weight of the 0 end, weight of the 1 end] per dimension).
 */
template<typename Float, std::size_t DIM, typename VT>
class RegionTrapezoidal {
public:
	using value_type = VT;
	static constexpr std::size_t dimensions = DIM;
	static constexpr std::size_t samples = 2;
	const Range<Float, DIM>& range() const { return _range; }

private:
	static constexpr std::size_t corners = std::size_t(1) << DIM;
	using Weights = std::array<std::array<Float,2>,DIM>;

	Range<Float, DIM> _range;
	std::array<value_type,corners> values; //Flat, the value of the corner m is values[m]

	void copy_corners(const multiarray<value_type,samples,DIM>& d) {
		for (std::size_t m = 0; m<corners; ++m) values[m] = d[corner_index(m)];
	}

	template<typename F>
	void evaluate_corners(const F& f) {
		const auto fr = this->f_in_range(f);
		for (std::size_t m = 0; m<corners; ++m) {
			std::array<double,DIM> p;
			for (std::size_t d = 0; d<DIM; ++d) p[d] = double((m >> d) & 1);
			values[m] = fr(p);
		}
	}

	template<typename F>
	constexpr auto f_in_range(const F& f) const { 
		return([&] (const std::array<double,DIM>& p) {
			std::array<Float,DIM> prange;
			for (std::size_t i = 0; i<DIM; ++i) 
				prange[i] = p[i]*(range().max(i) - range().min(i)) + range().min(i);
			return f(prange);
		}); 
	}

	constexpr Float volume_from(std::size_t start) const {
		Float v(1); 
		for (std::size_t i = start; i<DIM; ++i) v*=std::abs(range().max(i) - range().min(i));
		return v;	
	}

	//Corner m has the bit d of m as the index in dimension d
	static constexpr std::array<std::size_t,DIM> corner_index(std::size_t m) {
		std::array<std::size_t,DIM> idx{};
		for (std::size_t d = 0; d<DIM; ++d) idx[d] = (m >> d) & 1;
		return idx;
	}

	//Weights of the values at 0 and 1 for the integral of the interpolation between a and b 
	static constexpr std::array<Float,2> integral_weights(Float a, Float b) {
		Float w1 = (b*b - a*a)/Float(2);
		return {(b - a) - w1, w1};
	}

	static constexpr std::array<Float,2> mean_weights() { return {Float(0.5),Float(0.5)}; }

	static constexpr std::array<Float,corners> corner_weights(const Weights& w) {
		std::array<Float,corners> cw{};
		cw[0] = Float(1);
		for (std::size_t d = 0; d<DIM; ++d) {
			const std::size_t half = std::size_t(1) << d;
			for (std::size_t m = 0; m<half; ++m) {
				cw[m + half] = cw[m]*w[d][1];
				cw[m] *= w[d][0];
			}
		}
		return cw;
	}

	//Contracts one dimension at a time (the last one first), which keeps the dependency chains short 
	//and lets the compiler unroll everything for small dimensions
	value_type weighted_sum(const Weights& w) const {
		if constexpr (DIM == 0) return values[0];
		else {
			//The first contraction reads directly from the stored values
			constexpr std::size_t top = corners >> 1;
			std::array<value_type,top> v;
			for (std::size_t m = 0; m<top; ++m) 
				v[m] = values[m]*w[DIM-1][0] + values[m + top]*w[DIM-1][1];
			for (std::size_t d = DIM - 1; d-- > 0; ) {
				const std::size_t half = std::size_t(1) << d;
				for (std::size_t m = 0; m<half; ++m) 
					v[m] = v[m]*w[d][0] + v[m + half]*w[d][1];
			}
			return v[0];
		}
	}

public:
	RegionTrapezoidal(const std::array<Float, DIM>& range_min, 
			const std::array<Float, DIM>& range_max,
			const multiarray<value_type,samples,DIM>& d) : 
				_range(range_min, range_max) { copy_corners(d); }

	template<typename F>
	RegionTrapezoidal(const F& f, 
			const std::array<Float, DIM>& range_min, 
			const std::array<Float, DIM>& range_max) : 
				_range(range_min, range_max) {
		evaluate_corners(f);
	}
    
	RegionTrapezoidal(const Range<Float,DIM>& r, const multiarray<value_type,samples,DIM>& d) : 
				_range(r) { copy_corners(d); }

	template<typename F>
	RegionTrapezoidal(const F& f, const Range<Float,DIM>& r) : _range(r) {
		evaluate_corners(f);
	}

	value_type integral() const { 
		Weights w; w.fill(mean_weights());
		return range().volume()*weighted_sum(w); 
	}

	//If you do not reach the end it integrates in the other dimensions
	template<std::size_t DIMSUB>
	value_type approximation_at(const std::array<Float,DIMSUB>& pos) const {
		static_assert(DIM>=DIMSUB,"Cannot approximate with bigger number of dimensions than the region");
		const std::array<Float,DIMSUB> p = range().pos_in_range(pos);
		Weights w;
		for (std::size_t d = 0; d<DIM; ++d) 
			w[d] = (d<DIMSUB)?std::array<Float,2>{Float(1)-p[d],p[d]}:mean_weights();
		return weighted_sum(w)*volume_from(DIMSUB);
	}
	
	//If you do not reach the end it integrates in the other dimensions
	value_type approximation_at(Float pos) const {
		return approximation_at(std::array<Float,1>{pos});
	}

	//If you do not reach the end it integrates in the other dimensions
	template<std::size_t DIMSUB>
	value_type integral_subrange(const std::array<Float,DIMSUB>& a, 
						const std::array<Float,DIMSUB>& b) const {
		static_assert(DIM>=DIMSUB,"Cannot calculate the subrange integral for that many dimensions, as the region has less dimensions");
		const std::array<Float,DIMSUB> ua = range().pos_in_range(a), ub = range().pos_in_range(b);
		Weights w;
		for (std::size_t d = 0; d<DIM; ++d) 
			w[d] = (d<DIMSUB)?integral_weights(ua[d],ub[d]):mean_weights();
		return range().volume()*weighted_sum(w);
	}

	template<std::size_t DIMSUB>
	value_type integral_subrange(const Range<Float,DIMSUB>& subrange) const {
		return integral_subrange(subrange.min(),subrange.max());
	}
	
	value_type integral_subrange(Float a, Float b) const {
		return integral_subrange(std::array<Float,1>{a},std::array<Float,1>{b});
	}

private:
	//Inverse transform sampling of the line g(t) = g0*(1-t) + g1*t (non negative) in [a,b] 
	static Float sample_line(Float s, Float a, Float b, Float g0, Float g1) {
		const Float ga = g0*(Float(1)-a) + g1*a;
		const Float gb = g0*(Float(1)-b) + g1*b;
		const Float den = ga + std::sqrt((Float(1)-s)*ga*ga + s*gb*gb);
		const Float u = (den>Float(0))?(s*(b-a)*(ga+gb)/den):(s*(b-a)); //Uniform if the pdf is zero
		return std::min(b,std::max(a,a+u));
	}

public:
	/**
	 * Samples a point in the subrange [a,b] (in the first DIMSUB dimensions) with a pdf proportional 
	 * to the interpolation of the absolute value of the samples, integrating the remaining dimensions.
	 * The dimensions are sampled in order, each one from its marginal distribution given the previous ones.
	 * Only available when value_type is a floating point number.
	 */
	template<std::size_t DIMSUB, typename V = value_type, 
			typename = std::enable_if_t<std::is_floating_point_v<V>>>
	std::array<Float,DIMSUB> sample_subrange(const std::array<Float,DIMSUB>& s, 
						const std::array<Float,DIMSUB>& a, const std::array<Float,DIMSUB>& b) const {
		static_assert(DIM>=DIMSUB,"Cannot calculate the subrange sample for that many dimensions, as the region has less dimensions");
		const std::array<Float,DIMSUB> ua = range().pos_in_range(a), ub = range().pos_in_range(b);

		std::array<Float,corners> abs_data;
		for (std::size_t m = 0; m<corners; ++m) abs_data[m] = Float(std::abs(values[m]));

		Weights w;
		for (std::size_t d = 0; d<DIM; ++d) 
			w[d] = (d<DIMSUB)?integral_weights(ua[d],ub[d]):mean_weights();

		std::array<Float,DIMSUB> x;
		for (std::size_t d = 0; d<DIMSUB; ++d) {
			w[d] = {Float(1),Float(1)}; //Marginal in dimension d
			const std::array<Float,corners> cw = corner_weights(w);
			Float g0(0), g1(0);
			for (std::size_t m = 0; m<corners; ++m) 
				((m >> d) & 1 ? g1 : g0) += abs_data[m]*cw[m];
			x[d] = sample_line(s[d],ua[d],ub[d],g0,g1);
			w[d] = {Float(1)-x[d],x[d]}; //Fixed for the following dimensions
		}
		return range().pos_from_range(x);
	}

	template<std::size_t DIMSUB, typename V = value_type, 
			typename = std::enable_if_t<std::is_floating_point_v<V>>>
	std::array<Float,DIMSUB> sample_subrange(const std::array<Float,DIMSUB>& s, const Range<Float,DIMSUB>& subrange) const {
		return this->sample_subrange(s,subrange.min(),subrange.max());
	}
};

template<typename Float, typename F, std::size_t DIM>
auto region_trapezoidal(const F& f, 
			const std::array<Float, DIM>& range_min, 
			const std::array<Float, DIM>& range_max) {
	static_assert(DIM<=VILTRUM_MAX_DIMENSIONS_REGION,"Too many dimensions for region construction");
	return RegionTrapezoidal<Float,DIM,std::decay_t<decltype(f(range_min))>>(f,range_min, range_max);
}

template<typename Float, typename F, std::size_t DIM>
auto region_trapezoidal(const F& f, const Range<Float,DIM>& range) {
	static_assert(DIM<=VILTRUM_MAX_DIMENSIONS_REGION,"Too many dimensions for region construction");
	return RegionTrapezoidal<Float,DIM,std::decay_t<decltype(f(range.min()))>>(f,range);
}

}
