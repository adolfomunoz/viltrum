#pragma once

#include <array>
#include "../multiarray/multiarray.h"
#include "../multiarray/fill.h"
#include "../multiarray/fold.h"
#include "../multiarray/transform.h"
#include "../multiarray/io.h"
#include "../multiarray/split.h"
#include "../range.h"
#include "../control-variates/norm.h"
#include "../nested/nested.h"
#include <cmath>

namespace viltrum {


template<typename Float, typename Q, std::size_t DIM, typename VT>
class Region {
	Q quadrature;
	Range<Float, DIM> _range;
	
public:
	using value_type = VT;
	static constexpr std::size_t dimensions = DIM;
	using rule = Q;
	const Range<Float, DIM>& range() const {
		return _range;
	}
	
private:
	multiarray<value_type,Q::samples,DIM> data;
	
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
public:
	Region(const Q& q, 
			const std::array<Float, DIM>& range_min, 
			const std::array<Float, DIM>& range_max,
			multiarray<value_type,Q::samples,DIM>&& d) : 
				quadrature(q), _range(range_min, range_max), 
				data(std::forward<multiarray<value_type,Q::samples,DIM>>(d)) { }

	template<typename F>
	Region(const F& f, const Q& q, 
			const std::array<Float, DIM>& range_min, 
			const std::array<Float, DIM>& range_max) : 
				quadrature(q), _range(range_min, range_max) {
		data.fill(this->f_in_range(f));
	}
    
	Region(const Q& q, const Range<Float,DIM>& r, multiarray<value_type,Q::samples,DIM>&& d) : 
				quadrature(q), _range(r), 
				data(std::forward<multiarray<value_type,Q::samples,DIM>>(d)) { }

	template<typename F>
	Region(const F& f, const Q& q, const Range<Float,DIM>& r) : 
				quadrature(q), _range(r) {
		data.fill(this->f_in_range(f));
	}

    const Q& quadrature_rule() const { return quadrature; }
	
	value_type integral() const { return range().volume()*data.fold_all(quadrature); }
	
private:

	//Allegedly starts at 0
	template<typename MA, std::size_t DIMSUB>
	value_type app_at(const MA& ma, const std::array<Float,DIMSUB>& pos, 
			std::enable_if_t<MA::dimensions <= (DIM-DIMSUB+1)>* p = nullptr) const {
		return ma.fold([&] (const auto& v) 
				{ return quadrature.at(pos[DIMSUB - 1],v);}).fold_all(quadrature)*volume_from(DIMSUB);
	}
	
	template<typename MA, std::size_t DIMSUB>
	auto app_at(const MA& ma, const std::array<Float,DIMSUB>& pos, 
			std::enable_if_t<MA::dimensions >= ((DIM-DIMSUB)+2)>* p = nullptr) const {
		return app_at(ma.fold([&] (const auto& v) 
				{ return quadrature.at(pos[DIMSUB - MA::dimensions],v); }),pos);
	}

public:
	//If you do not reach the end it integrates in the other dimensions
	template<std::size_t DIMSUB>
	value_type approximation_at(const std::array<Float,DIMSUB>& pos) const {
		static_assert(DIM>=DIMSUB,"Cannot approximate with bigger number of dimensions than the region");
		return app_at(data,range().pos_in_range(pos));
	}
	
	//If you do not reach the end it integrates in the other dimensions
	value_type approximation_at(Float pos) const {
		return approximation_at(std::array<Float,1>{pos});
	}

private:
	//Starts from 0
	template<typename MA, std::size_t DIMSUB>
	value_type sub_first(const MA& ma,  
			const std::array<Float,DIMSUB>& a, const std::array<Float,DIMSUB>& b, 
			std::enable_if_t<MA::dimensions <= (DIM-DIMSUB+1)>* p = nullptr) const {
		return ma.fold([&] (const auto& v) 
			{ return quadrature.subrange(a[DIMSUB - 1], b[DIMSUB - 1],v); }).fold_all(quadrature);
	}
	
	//Allegedly starts at 0
	template<typename MA, std::size_t DIMSUB>
	auto sub_first(const MA& ma, 
			const std::array<Float,DIMSUB>& a, const std::array<Float,DIMSUB>& b, 
			std::enable_if_t<MA::dimensions >= ((DIM-DIMSUB)+2)>* p = nullptr) const {
		return sub_first(ma.fold([&] (const auto& v) 
			{ return quadrature.subrange(a[DIMSUB - MA::dimensions], b[DIMSUB- MA::dimensions],v);}),a,b);
	}
public:
	template<std::size_t DIMSUB>
	value_type integral_subrange_first(const std::array<Float,DIMSUB>& a, 
						const std::array<Float,DIMSUB>& b) const {
		static_assert(DIM>=DIMSUB,"Cannot calculate the subrange integral for that many dimensions, as the region has less dimensions");
		return range().volume()*sub_first(data,range().pos_in_range(a),range().pos_in_range(b));
	}

private:
	template<typename MA, std::size_t DIMSUB>
	value_type sub_last(const MA& ma, 
			const std::array<Float,DIMSUB>& a, const std::array<Float,DIMSUB>& b) const {
        if constexpr (MA::dimensions > DIMSUB)
            return sub_last(ma.fold(quadrature,DIMSUB),a,b);
        else if constexpr (MA::dimensions > 1)
            return sub_last(ma.fold([&] (const auto& v) { return quadrature.subrange(a[MA::dimensions-1],b[MA::dimensions-1],v); },MA::dimensions-1),a,b);
        else
            return ma.fold([&] (const auto& v) { return quadrature.subrange(a[0],b[0],v); }).value();
	}
public:
	template<std::size_t DIMSUB>
	value_type integral_subrange_last(const std::array<Float,DIMSUB>& a, 
						const std::array<Float,DIMSUB>& b) const {
		static_assert(DIM>=DIMSUB,"Cannot calculate the subrange integral for that many dimensions, as the region has less dimensions");
		return range().volume()*sub_last(data,range().pos_in_range(a),range().pos_in_range(b));
    }

   
	template<std::size_t DIMSUB>
	value_type integral_subrange(const std::array<Float,DIMSUB>& a, 
						const std::array<Float,DIMSUB>& b) const {
        return integral_subrange_last(a,b);
	}

    template<std::size_t DIMSUB>
    value_type integral_subrange(const Range<Float,DIMSUB>& subrange) const {
        return integral_subrange(subrange.min(),subrange.max());
    }
	
	value_type integral_subrange(Float a, Float b) const {
		return integral_subrange(std::array<Float,1>{a},std::array<Float,1>{b});
	}

private:
/*
	template<typename MA, std::size_t DIMSUB, typename Norm = NormDefault>
	Float sample_sub(const MA& ma, 
			const Float& s, const std::array<Float,DIMSUB>& a, const std::array<Float,DIMSUB>& b,
			const Norm& norm = Norm()) const {
        if constexpr (MA::dimensions > DIMSUB)
            return sample_sub(ma.fold(
				[&] (const auto& v) { return std::max(Float(1.e-20),quadrature.pdf_integral_subrange(a[DIMSUB],b[DIMSUB],v,norm)); },DIMSUB),
				s,a,b);
        else if constexpr (MA::dimensions > 1)
			return sample_sub(ma.fold([&] (const auto& v) {
				return std::max(Float(1.e-20),quadrature.pdf_integral_subrange(a[0],b[0],v));
			},0),s,pop(a),pop(b));
		else
            return ma.fold([&] (const auto& v) { 
				return quadrature.sample(s,a[0],b[0],v,norm); },0).value();
	}

	template<typename MA, std::size_t DIMSUB, typename Norm = NormDefault>
	std::array<Float,DIMSUB> sample_subrange_i(const MA& ma, 
			const std::array<Float,DIMSUB>& sol_i,
			const std::array<Float,DIMSUB>& s, const std::array<Float,DIMSUB>& a, 
			const std::array<Float,DIMSUB>& b,
			const Norm& norm = Norm()) const {
        if constexpr (MA::dimensions > DIMSUB)
            return sample_subrange_i(ma.fold(
				[&] (const auto& v) { return quadrature.pdf_integral_subrange(a[DIMSUB],b[DIMSUB],v,norm); },DIMSUB),
				sol_i,s,a,b,norm);
		else { 
			std::array<Float,DIMSUB> sol = sol_i;
			sol[MA::dimensions-1] = sample_sub(ma,s[MA::dimensions-1],a,b,norm); 
			if constexpr (MA::dimensions == 1)
				return sol;
			else 
				return sample_subrange_i(
					ma.fold([&] (const auto& v) {
						return quadrature.pdf(sol[MA::dimensions-1],v,a[MA::dimensions-1],b[MA::dimensions-1],norm);
//						return std::max(Float(1.e-20),norm(quadrature.at(sol[MA::dimensions-1],v)));
//						return quadrature.pdf(sol[MA::dimensions-1],v,a[MA::dimensions-1],b[MA::dimensions-1]);	
					},MA::dimensions-1),sol,s,a,b,norm);
		}	
	} 
*/

	template<typename MA, std::size_t DIMSUB, typename Norm = NormDefault>
	Float sample_marginal(const MA& ma, 
			const std::array<Float,DIMSUB>& s, const std::array<Float,DIMSUB>& a, 
			const std::array<Float,DIMSUB>& b,
			const Norm& norm = Norm()) const {

		if constexpr (MA::dimensions == 1) {
			return ma.fold([&] (const auto& v) { 
				return quadrature.sample(s[0],a[0],b[0],v,norm); }).value();
		} else if constexpr (DIMSUB < MA::dimensions) {
			return sample_marginal(ma.fold([&] (const auto& v) {
				return quadrature.pdf_integral_subrange(Float(0),Float(1),v,norm);
				}, MA::dimensions-1),s,a,b,norm);
		} else {
			return sample_marginal(ma.fold([&] (const auto& v) {
				return quadrature.pdf_integral_subrange(a[1],b[1],v,norm);
				}, MA::dimensions-1),s,a,b,norm);
		} 
	}

	template<typename MA, std::size_t DIMSUB, typename Norm = NormDefault>
	std::array<Float,DIMSUB> sample_subrange_i(const MA& ma,
			const std::array<Float,DIMSUB>& s, const std::array<Float,DIMSUB>& a, 
			const std::array<Float,DIMSUB>& b,
			const Norm& norm = Norm()) const {
		
		Float sol0 = sample_marginal(ma,s,a,b,norm);
		if constexpr (DIMSUB==1) return std::array<Float,1>{sol0};
		else return sol0|sample_subrange_i(
			ma.fold([&] (const auto& v) {
				return quadrature.pdf_unnormalized(sol0,v,a[0],b[0],norm);
			}),pop(s),pop(a),pop(b),norm); 	
	} 

public:
	template<std::size_t DIMSUB,typename Norm = NormDefault>
	std::array<Float,DIMSUB> sample_subrange(const std::array<Float,DIMSUB>& s, const std::array<Float,DIMSUB>& a, 
						const std::array<Float,DIMSUB>& b, const Norm& norm = Norm()) const {
		static_assert(DIM>=DIMSUB,"Cannot calculate the subrange sample for that many dimensions, as the region has less dimensions");
		return range().pos_from_range(
			sample_subrange_i(data,s,range().pos_in_range(a),range().pos_in_range(b),norm));
	} 

	template<std::size_t DIMSUB,typename Norm = NormDefault>
	std::array<Float,DIMSUB> sample_subrange(const std::array<Float,DIMSUB>& s, const Range<Float,DIMSUB>& range, const Norm& norm = Norm()) const {
		return this->sample_subrange(s,range.min(),range.max(),norm);
	} 

private:
	template<typename MA, std::size_t DIMSUB, typename Norm = NormDefault>
	Float pdf_sub(const MA& ma, 
			const std::array<Float,DIMSUB>& a, const std::array<Float,DIMSUB>& b,
			const Norm& norm = Norm()) const {
		if constexpr (MA::dimensions == 1) return ma.fold([&] (const auto& v) {
			return quadrature.pdf_integral_subrange(a[0],b[0],v,norm);
			}).value();
		else if constexpr (DIMSUB<MA::dimensions) return pdf_sub(ma.fold([&] (const auto& v) {
				return quadrature.pdf_integral_subrange(Float(0),Float(1),v,norm);
			},MA::dimensions-1),a,b,norm);
		else return pdf_sub(ma.fold([&] (const auto& v) {
				return quadrature.pdf_integral_subrange(a[MA::dimensions-1] ,b[MA::dimensions-1] ,v,norm);
			},MA::dimensions-1),a,b,norm);
	}

public:
	template<std::size_t DIMSUB, typename Norm = NormDefault>
	Float pdf_integral_subrange(const std::array<Float,DIMSUB>& a, 
						const std::array<Float,DIMSUB>& b, const Norm& norm = Norm()) const {
		static_assert(DIM>=DIMSUB,"Cannot calculate the subrange integral for that many dimensions, as the region has less dimensions");
		return range().volume()*pdf_sub(data,range().pos_in_range(a),range().pos_in_range(b),norm);
	}

	template<std::size_t DIMSUB,typename Norm = NormDefault>
	Float pdf_integral_subrange(const Range<Float,DIMSUB>& range, const Norm& norm = Norm()) const {
		return this->pdf_integral_subrange(range.min(),range.max(),norm);
	} 



private:
	template<typename MA, std::size_t DIMSUB, typename Norm = NormDefault>
	Float pdf_at(const MA& ma, const std::array<Float,DIMSUB>& pos, const std::array<Float,DIMSUB>& a, 
						const std::array<Float,DIMSUB>& b, const Norm& norm = Norm()) const {
		if constexpr (MA::dimensions == 1) return ma.fold([&] (const auto& v) {
			return quadrature.pdf_unnormalized(pos[0],v,a[0],b[0],norm);
			}).value();
		else if constexpr (DIMSUB<MA::dimensions) return pdf_at(ma.fold([&] (const auto& v) {
				return quadrature.pdf_integral_subrange(Float(0),Float(1),v,norm);
			},MA::dimensions-1),pos,a,b,norm);
		else return pdf_at(ma.fold([&] (const auto& v) {
				return quadrature.pdf_unnormalized(
					pos[MA::dimensions-1],v,a[MA::dimensions-1],b[MA::dimensions-1],norm);
			},MA::dimensions-1),pos,a,b,norm);

	}

	

public:

#if 0
	template<std::size_t DIMSUB,typename Norm = NormDefault>
	Float pdf_subrange(const std::array<Float,DIMSUB>& pos, const std::array<Float,DIMSUB>& a, 
						const std::array<Float,DIMSUB>& b,const Norm& norm = Norm()) const {
		static_assert(DIM>=DIMSUB,"Cannot calculate the subrange pdf for that many dimensions, as the region has less dimensions");
		Float den = this->pdf_integral_subrange(a,b,norm);
		if (den<1.e-3) return Float(1)/viltrum::range(a,b).volume();
		else return std::max(Float(0),norm(this->approximation_at(pos))/den);
	}
#else
    template<std::size_t DIMSUB,typename Norm = NormDefault>
	Float pdf_subrange(const std::array<Float,DIMSUB>& pos, const std::array<Float,DIMSUB>& a, 
						const std::array<Float,DIMSUB>& b,const Norm& norm = Norm()) const {
		static_assert(DIM>=DIMSUB,"Cannot calculate the subrange pdf for that many dimensions, as the region has less dimensions");
		return pdf_at(data,range().pos_in_range(pos),range().pos_in_range(a),range().pos_in_range(b),norm)/(pdf_integral_subrange(a,b,norm));
	}
#endif

	template<std::size_t DIMSUB,typename Norm = NormDefault>
	Float pdf_subrange(const std::array<Float,DIMSUB>& pos, const Range<Float,DIMSUB>& range, const Norm& norm = Norm()) const {
		return this->pdf_subrange(pos,range.min(),range.max(),norm);
	} 
	
	template<typename F>
	std::vector<Region<Float,Q,DIM,value_type>> split(const F& f, std::size_t dimension = 0, std::size_t parts = 2) const {
		assert(dimension < DIM);
		auto subdatas = data.split(f_in_range(f),dimension,parts);
		std::vector<Region<Float,Q,DIM,value_type>> sol; sol.reserve(parts);
		std::array<Float,DIM> range_midmin = range().min();
		std::array<Float,DIM> range_midmax = range().max();
		Float d = (range().max(dimension) - range().min(dimension))/(Float(parts));
		for (std::size_t i = 0; i<parts; ++i) {
			range_midmax[dimension] = (i<(parts-1))?(range().min(dimension)+d*(i+1)):range().max(dimension);
			sol.emplace_back(quadrature,range_midmin,range_midmax,std::move(subdatas[i]));
			range_midmin[dimension] = range_midmax[dimension];
		}
		return sol;
	}
	
	template<typename F>
	std::vector<Region<Float,Q,DIM,value_type>> split_all(const F& f, std::size_t parts = 2) const {
		std::vector<Region<Float,Q,DIM,value_type>> sol, accumulated;
		sol.reserve(std::size_t(std::pow(parts,DIM))); 
		sol.push_back(*this);
		for (std::size_t d = 0; d < DIM; ++d) {
			accumulated.clear();
			accumulated.reserve(sol.capacity());
			for (auto r : sol) 
				for (auto subregion : r.split(f,d,parts))
					accumulated.emplace_back(subregion);
			
			sol.swap(accumulated);
		}
		return sol;
	}
	

	template<typename QDT = Q, typename = typename std::enable_if<is_nested<QDT>::value>::type >
	value_type error(std::size_t dim) const {
		assert(dim < DIM);
		return range().volume()*
			data.fold([&] (const auto& v) { return quadrature.error(v); },dim)
			    .fold_all(quadrature);
	}

	template<typename ErrorMetric, typename QDT = Q, typename = typename std::enable_if<is_nested<QDT>::value>::type >
	auto error(std::size_t dim,const ErrorMetric& error_metric) const {
		assert(dim < DIM);
		return range().volume()*
			data.fold([&] (const auto& v) { return quadrature.error(v,error_metric); },dim)
			    .fold_all(quadrature);
	}

	
	/*
	 * ErrorMetric :: (value_type,value_type) -> arithmetic
	 */


	template<typename ErrorMetric, typename QDT = Q, typename = typename std::enable_if<is_nested<QDT>::value>::type >
	auto max_error_dimension(const ErrorMetric& error_metric) const {
		Float max_err = 0; std::size_t max_dim = 0; Float err;
		for (std::size_t d = 0; d<DIM; ++d) {
			err = this->error(d,error_metric);
			if (err>max_err) {
				max_err = err; max_dim = d;
			}
		}
		return std::tuple(max_err,max_dim);
	}

	template<typename QDT = Q, typename = typename std::enable_if<is_nested<QDT>::value && std::is_arithmetic<value_type>::value>::type >
	std::tuple<Float,std::size_t> max_error_dimension() const {
		return max_error_dimension(
			[] (const value_type& v1,const value_type& v2) { return Float(std::abs(v2-v1)); }
		);
	}

	template<typename QDT = Q, typename = typename std::enable_if<is_nested<QDT>::value>::type >
	value_type error() const {
		return range().volume()*(data.fold_all(quadrature) - 
					   data.fold_all([&] (const auto& v) { return quadrature.low(v); }));
	}		
};

template<typename Float, typename F, typename Q, std::size_t DIM>
auto region(const F& f, const Q& q, 
			const std::array<Float, DIM>& range_min, 
			const std::array<Float, DIM>& range_max) {
	return Region<Float,Q,DIM,std::decay_t<decltype(f(range_min))>>(f,q,range_min, range_max);
}

template<typename Float, typename F, typename Q, std::size_t DIM>
auto region(const F& f, const Q& q, const Range<Float,DIM>& range) {
	return Region<Float,Q,DIM,std::decay_t<decltype(f(range.min()))>>(f,q,range);
}

template<typename R, typename E>
class ExtendedRegion : public R {
	E e;
public:
	ExtendedRegion(R&& r, E&& e) : R(std::forward<R>(r)), e(std::forward<E>(e)) { }
	ExtendedRegion(const R& r, E&& e) : R(r), e(std::forward<E>(e)) { }
	ExtendedRegion(R&& r, const E& e) : R(std::forward<R>(r)), e(e) { }
	ExtendedRegion(const R& r, const E& e) : R(r), e(e) { }
	
	const E& extra() const { return e; }
};

template<typename Float, std::size_t DIM, std::size_t DIMBINS, typename R>
auto pixels_in_region(const R& r, const std::array<std::size_t,DIMBINS>& bin_resolution, const Range<Float,DIM>& range) {
	std::array<std::size_t,DIMBINS> start_bin, end_bin;
    for (std::size_t i = 0; i<DIMBINS;++i) {
        start_bin[i] = std::max(std::size_t(0),std::size_t(Float(bin_resolution[i])*(r.range().min(i) - range.min(i))/(range.max(i) - range.min(i))));
        end_bin[i]   = std::max(start_bin[i]+1,std::min(bin_resolution[i],std::size_t(0.99f + 
                    (Float(bin_resolution[i])*(r.range().max(i) - range.min(i))/(range.max(i) - range.min(i))))));
    }
	return multidimensional_range(start_bin,end_bin);
}

}
