#pragma once
#include <random>
#include "../range-infinite.h"


namespace viltrum {


//I am implementing it so it is predictable once we initialize the RNG.
template<typename RNG, typename Number>
class RandomSequenceRNG {
    RNG rng;
    RangeInfinite<Number> range;
    
public:
    RandomSequenceRNG(RNG&& r, const RangeInfinite<Number>& ra = RangeInfinite<Number>()) : 
        rng(std::forward<RNG>(r)), range(ra) {}

    class const_iterator {
    private:
        RNG rng;
        const RangeInfinite<Number>& range;
        std::size_t i;
        std::uniform_real_distribution<Number> dis;
        Number n;
        const_iterator(const RNG& r, const RangeInfinite<Number>& ra) :
            rng(r), range(ra), i(0), n(std::uniform_real_distribution<Number>(range.min(i),range.max(i))(rng)) {} 
        friend class RandomSequenceRNG<RNG,Number>;
    public:
        const Number& operator*() const { return n; }
        const_iterator& operator++()    {  ++i; n = std::uniform_real_distribution<Number>(range.min(i),range.max(i))(rng); return (*this); }
		bool operator==(const const_iterator& that) const { return false; } //Infinite list
		bool operator!=(const const_iterator& that) const { return true;  } //Infinite list
    }; 

    const_iterator begin() const { return const_iterator(rng,range); }
    const_iterator end() const { return const_iterator(rng,range); }
};

template<typename Number,typename RNG = std::mt19937>
auto random_sequence_rng(const RangeInfinite<Number>& r, std::size_t seed = std::random_device()()) {
    return RandomSequenceRNG<RNG,Number>(RNG(seed),r);
} 

}
