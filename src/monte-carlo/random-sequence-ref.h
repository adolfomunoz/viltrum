#pragma once
#include <random>
#include "../range-infinite.h"


namespace viltrum {


//I am implementing it so it is predictable once we initialize the RNG.
template<typename RNG, typename Number>
class RandomSequenceRef {
    RNG& rng;
    RangeInfinite<Number> range;
    
public:
    RandomSequenceRef(RNG& r, const RangeInfinite<Number>& ra = RangeInfinite<Number>()) : 
        rng(r), range(ra) {}

    class const_iterator {
    private:
        RNG& rng;
        const RangeInfinite<Number>& range;
        std::size_t i;
        std::uniform_real_distribution<Number> dis;
        Number n;
        const_iterator(RNG& r, const RangeInfinite<Number>& ra) :
            rng(r), range(ra), i(0), n(std::uniform_real_distribution<Number>(range.min(i),range.max(i))(rng)) {} 
        friend class RandomSequenceRef<RNG,Number>;
    public:
        const Number& operator*() const { return n; }
        const_iterator& operator++()    {  ++i; n = std::uniform_real_distribution<Number>(range.min(i),range.max(i))(rng); return (*this); }
		bool operator==(const const_iterator& that) const { return false; } //Infinite list
		bool operator!=(const const_iterator& that) const { return true;  } //Infinite list
    }; 

    const_iterator begin() const { return const_iterator(rng,range); }
    const_iterator end() const { return const_iterator(rng,range); }
};

template<typename Number,typename RNG>
auto random_sequence_ref(const RangeInfinite<Number>& r, RNG& rng) {
    return RandomSequenceRef<RNG,Number>(rng,r);
} 

}
