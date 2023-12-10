#pragma once
#include <array>
#include "../newton-cotes/newton-cotes.h"

namespace viltrum {

template<typename H, typename L>
class Nested : public H {
	L l;
public:	
	static constexpr std::size_t samples = H::samples;
	static_assert(((H::samples-1)%(L::samples-1))==0, "Quadrature rules are not compatible for a nested rule");
	
	Nested(const H& h, const L& l) : H(h), l(l) { }
	
	template<typename T>
	auto low(const std::array<T,samples>& p) const -> T {
		std::array<T,L::samples> plow;
		for (std::size_t i = 0; i < L::samples; ++i) {
			plow[i] = p[(i*(H::samples-1)/(L::samples-1))];
		}
		return l(plow);
	}
	
	template<typename T>
	auto error(const std::array<T,samples>& p) const -> T {
		return (*this)(p) - low(p);
	}

	template<typename T, typename ErrorMetric>
	auto error(const std::array<T,samples>& p, const ErrorMetric& error_metric) const {
		return error_metric((*this)(p),low(p));
	}
};

template<typename T>
struct is_nested : std::integral_constant<bool,false> {};
template<typename H,typename L>
struct is_nested<Nested<H,L>> : std::integral_constant<bool,true> {};


template<typename H, typename L>
auto nested(const H& h, const L& l) {
	return Nested<H,L>(h,l);
}

}
