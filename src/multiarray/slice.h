#pragma once
#include "multiarray-crtp.h"
#include <cassert>

namespace viltrum {
    
//Class that fixes one dimension to a specific index and returns the rest
//static_assert when base dimension is 1
template<typename Base>
class multiarray_slice_const : public multiarray_const<multiarray_slice_const<Base>> {
	const Base& base;
	std::size_t dimension;
	std::size_t index;
public:
	static constexpr std::size_t size = Base::size;
	static constexpr std::size_t dimensions = Base::dimensions - 1;
	using value_type = typename Base::value_type;
	using index_type = std::array<std::size_t,dimensions>;


	multiarray_slice_const(const Base& base, std::size_t dimension = 0, std::size_t index = 0) noexcept :base(base),dimension(dimension),index(index) {}
	const value_type& operator[](const index_type& indices) const noexcept { 
		std::array<std::size_t,Base::dimensions> slice_indices; 
		for (std::size_t d = 0; d < dimension; ++d)
			slice_indices[d] = indices[d];
		slice_indices[dimension]=index;
		for (std::size_t d = dimension; d < dimensions; ++d)
			slice_indices[d + 1] = indices[d];
		return base[slice_indices]; 
	}
};

//Class that fixes one dimension to a specific index and returns the rest
template<typename Base>
class multiarray_slice : public multiarray_mutable<multiarray_slice<Base>> {
	Base& base;
	std::size_t  dimension;
	std::size_t index;
public:
	static constexpr std::size_t size = Base::size;
	static constexpr std::size_t  dimensions = Base::dimensions - 1;
	using value_type = typename Base::value_type;
	using index_type = std::array<std::size_t,dimensions>;


	multiarray_slice(Base& base,std::size_t  dimension = 0, std::size_t index = 0) noexcept :base(base),dimension(dimension),index(index) {}
	const value_type& operator[](const index_type& indices) const noexcept { 
		std::array<std::size_t,Base::dimensions> slice_indices; 
		for (std::size_t d = 0; d < dimension; ++d)
			slice_indices[d] = indices[d];
		slice_indices[dimension]=index;
		for (std::size_t d = dimension; d < dimensions; ++d)
			slice_indices[d + 1] = indices[d];
		return base[slice_indices]; 
	}
	value_type& operator[](const index_type& indices) noexcept { 
		std::array<std::size_t,Base::dimensions> slice_indices; 
		for (std::size_t d = 0; d < dimension; ++d)
			slice_indices[d] = indices[d];
		slice_indices[dimension]=index;
		for (std::size_t d = dimension; d < dimensions; ++d)
			slice_indices[d + 1] = indices[d];
		return base[slice_indices]; 
	}
	
	ASSIGNMENT(multiarray_slice);
};

namespace detail {
template<typename MA>
auto slice(const MA& ma, std::size_t  dimension, std::size_t index) noexcept {
	static_assert(MA::dimensions > 1, "Cannot slice a dimension 1 multiarray");
	assert( (dimension >= 0) && (dimension<MA::dimensions) && (index < MA::size) );
	return multiarray_slice_const<MA>(ma,dimension,index);
}

template<typename MA>
auto slice(MA& ma, std::size_t  dimension, std::size_t index) noexcept {
	static_assert(MA::dimensions > 1, "Cannot slice a dimension 1 multiarray");
	assert( (dimension >= 0) && (dimension<MA::dimensions) && (index < MA::size) );	return multiarray_slice<MA>(ma,dimension,index);
}
}

}