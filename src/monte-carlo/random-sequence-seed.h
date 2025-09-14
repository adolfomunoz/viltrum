#pragma once
#include <random>
#include "../range-infinite.h"


namespace viltrum {


//I am implementing it so it is predictable once we initialize the RNG.
template<typename RNG, typename Number>
class RandomSequenceSeed {
    std::size_t seed;
    RangeInfinite<Number> range;
    
public:
    RandomSequenceSeed(std::size_t s, const RangeInfinite<Number>& ra = RangeInfinite<Number>()) : 
        seed(s), range(ra) {}

    class const_iterator {
    private:
        RNG rng;
        const RangeInfinite<Number>& range;
        std::size_t i;
        std::uniform_real_distribution<Number> dis;
        Number n;
        const_iterator(std::size_t seed, const RangeInfinite<Number>& ra) :
            rng(seed), range(ra), i(0), n(std::uniform_real_distribution<Number>(range.min(i),range.max(i))(rng)) {} 
        friend class RandomSequenceSeed<RNG,Number>;
    public:
        const Number& operator*() const { return n; }
        const_iterator& operator++()    {  ++i; n = std::uniform_real_distribution<Number>(range.min(i),range.max(i))(rng); return (*this); }
		bool operator==(const const_iterator& that) const { return false; } //Infinite list
		bool operator!=(const const_iterator& that) const { return true;  } //Infinite list
    }; 

    const_iterator begin() const { return const_iterator(seed,range); }
    const_iterator end() const { return const_iterator(0,range); }
};

template<typename Number,typename RNG = std::mt19937>
auto random_sequence_seed(const RangeInfinite<Number>& r, std::size_t seed = std::random_device()()) {
    return RandomSequenceSeed<RNG,Number>(seed,r);
} 


}
