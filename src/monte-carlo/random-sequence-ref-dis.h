#pragma once
#include <random>
#include "../range-infinite.h"


namespace viltrum {


//I am implementing it so it is predictable once we initialize the RNG.
template<typename RNG, typename Number>
class RandomSequenceRefDis {
    RNG& rng;
    mutable std::uniform_real_distribution<Number> dis;
    RangeInfinite<Number> range;
    
public:
    RandomSequenceRefDis(RNG& r, const RangeInfinite<Number>& ra = RangeInfinite<Number>()) : 
        rng(r), dis(Number(0),Number(1)),range(ra) {}

    class const_iterator {
    private:
        RNG& rng;
        const RangeInfinite<Number>& range;
        std::size_t i;
        std::uniform_real_distribution<Number>& dis;
        Number n;
        const_iterator(RNG& r, std::uniform_real_distribution<Number>& d, const RangeInfinite<Number>& ra) :
            rng(r), range(ra), i(0), dis(d), n(dis(rng)*(range.max(i)-range.min(i))+range.min(i)) {} 
        friend class RandomSequenceRefDis<RNG,Number>;
    public:
        const Number& operator*() const { return n; }
        const_iterator& operator++()    {  ++i; n = dis(rng)*(range.max(i)-range.min(i))+range.min(i); return (*this); }
		bool operator==(const const_iterator& that) const { return false; } //Infinite list
		bool operator!=(const const_iterator& that) const { return true;  } //Infinite list
    }; 

    const_iterator begin() const { return const_iterator(rng,dis,range); }
    const_iterator end() const { return const_iterator(rng,dis,range); }
};

template<typename Number,typename RNG>
auto random_sequence_ref_dis(const RangeInfinite<Number>& r, RNG& rng) {
    return RandomSequenceRefDis<RNG,Number>(rng,r);
} 

}
