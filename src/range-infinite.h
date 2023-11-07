#pragma once

#include <vector>
#include "range.h"

#if (__cplusplus < 201703L)
namespace std {
    template< class T >
    inline constexpr bool is_floating_point_v = std::is_floating_point<T>::value;
}
#endif

namespace viltrum {

template<typename T>
class RangeInfinite : public std::array<std::vector<T>,2> {
	T _volume;
public:
	RangeInfinite(const std::vector<T>& a = std::vector<T>(), const std::vector<T>& b = std::vector<T>()) :
		std::array<std::vector<T>,2>{a,b}  {
		_volume = T(1);
		for (std::size_t i = 0; i<std::max(a.size(),b.size()); ++i) _volume*=(max(i)-min(i));
	}

    using value_type = T;
//    static constexpr std::size_t size = DIM;
    
	const std::vector<T>& min() const { return (*this)[0]; }
	T min(std::size_t i) const {  
        return (i<this->min().size())?(this->min()[i]):T(0);
    }
	const std::vector<T>& max() const { return (*this)[1]; }
	T max(std::size_t i) const {  
        return (i<this->max().size())?(this->max()[i]):T(1);
    }
	T volume() const { return _volume; }
    
    template<typename C>
    bool is_inside(const C& x) const {
        bool is = true;
        for (std::size_t i = 0; (i<x.size()) && is; ++i)
           is = ( (x[i]>=min(i)) && (x[i]<=max(i)) );
        return is;
    }
	
	template<std::size_t DIMSUB>
	std::array<T,DIMSUB> pos_in_range(const std::array<T,DIMSUB>& pos) const {
		std::array<T,DIMSUB> prange;
		for (std::size_t i = 0; i<DIMSUB; ++i) 
			prange[i] = (pos[i]-min(i))/(max(i) - min(i));
		return prange;
	}
	
	RangeInfinite<T> subrange_dimension(std::size_t dim, T a, T b) const {
		std::vector<T> new_a = min();
        if (dim>=new_a.size()) new_a.resize(dim+1,T(0)); 
        new_a[dim]=a;
		std::vector<T> new_b = max(); 
        if (dim>=new_b.size()) new_b.resize(dim+1,T(1));
        new_b[dim]=b;
		return RangeInfinite<T>(new_a, new_b);
	}

    template<std::size_t DIMSUB>
    Range<T,DIMSUB> intersection(const Range<T,DIMSUB>& that) const {
		std::array<T,DIMSUB> a, b; 
        for (std::size_t d = 0; d<DIMSUB; ++d) {
            a[d] = std::max(this->min(d),that.min(d));
            b[d] = std::max(a[d],std::min(this->max(d),that.max(d))); //std::max watches out for empty ranges
        }
        return Range<T,DIMSUB>(a,b);
    }

    template<std::size_t DIMSUB>
    RangeInfinite<T> intersection_large(const Range<T,DIMSUB>& that) const {
        std::vector<T> a = this->min(), b = this->max();
        if (a.size()<DIMSUB) a.resize(DIMSUB);
        if (b.size()<DIMSUB) b.resize(DIMSUB);
        for (std::size_t d = 0; d<DIMSUB; ++d) {
            a[d] = std::max(this->min(d),that.min(d));
            b[d] = std::max(a[d],std::min(this->max(d),that.max(d))); //std::max watches out for empty ranges
        }
        return RangeInfinite<T>(a,b);
    }
};


template<typename T>
RangeInfinite<T> range_infinite(const std::vector<T>& a, const std::vector<T>& b) {
    return RangeInfinite<T>(a,b);
}

template<typename T, std::size_t N1>
auto operator|(const Range<T,N1>& r1,const RangeInfinite<T>& r2) noexcept -> RangeInfinite<T> {
    std::vector<T> a(r1.min().size()+r2.min().size());
    std::copy(r1.min().begin(),r1.min().end(),a.begin());
    std::copy(r2.min().begin(),r2.min().end(),a.begin()+N1);
    std::vector<T> b(r1.max().size()+r2.max().size());
    std::copy(r1.max().begin(),r1.max().end(),b.begin());
    std::copy(r2.max().begin(),r2.max().end(),b.begin()+N1);    
    return range_infinite(a,b);
}

template<typename T = float>
RangeInfinite<T> range_primary_infinite() {
    return RangeInfinite<T>();
}

}
