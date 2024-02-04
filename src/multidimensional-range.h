#pragma once
#include <array>

namespace viltrum {

template<std::size_t DIM>
class MultidimensionalRange {
    std::array<std::size_t,DIM> min_;
    std::array<std::size_t,DIM> max_;
    friend class const_iterator;

public:
    using value_type = std::array<std::size_t,DIM>;

    constexpr std::size_t min(std::size_t i) const { return min_[i]; }
    constexpr std::size_t max(std::size_t i) const { return max_[i]; }

    constexpr std::size_t size() const {
		std::size_t n = 1;
		for (std::size_t i = 0; i<DIM; ++i) n*=(max(i)-min(i));
		return n;
	}

    operator std::size_t() const { return size(); }

    class const_iterator {
        const MultidimensionalRange<DIM>& range;
        std::array<std::size_t,DIM> a; bool ended;
    public:
        const_iterator(const MultidimensionalRange<DIM>& range, const std::array<std::size_t,DIM>& a) : range(range), a(a),ended(false) {}
        const_iterator(const MultidimensionalRange<DIM>& range) : range(range), ended(true) {}
        const std::array<std::size_t,DIM>& operator*() const { return a; }
        const_iterator& operator++()    {
            for (std::size_t d = 0; d<DIM; ++d) {
                ++a[d];
                if (a[d]>=range.max(d)) a[d] = range.min(d);
                else return (*this);   
            }
            ended = true; return (*this);
        }
        const_iterator operator++(int) {
            const_iterator old = (*this);
            ++(*this);
            return old;
        }
        bool operator==(const const_iterator& that) { return (this->ended == that.ended); }
        bool operator!=(const const_iterator& that) { return (this->ended != that.ended); }
    };

    const_iterator begin() const { return const_iterator(*this, min_); }
    const_iterator end()   const { return const_iterator(*this); }

    MultidimensionalRange(const std::array<std::size_t,DIM>& a, const std::array<std::size_t,DIM>& b) :
        min_(a), max_(b)  { }
};

template<std::size_t DIM>
MultidimensionalRange<DIM> multidimensional_range(const std::array<std::size_t,DIM>& a, const std::array<std::size_t,DIM>& b) {
    return MultidimensionalRange<DIM>(a,b);
}


template<std::size_t DIM>
MultidimensionalRange<DIM> multidimensional_range(const std::array<std::size_t,DIM>& b) {
    std::array<std::size_t,DIM> a; for (std::size_t d = 0; d<DIM; ++d) a[d]=0;
    return multidimensional_range(a,b);
}




}



