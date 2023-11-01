#pragma once
#include <random>

namespace viltrum {

//I am implementing it so it is predictable once we initialize the RNG.
template<typename RNG, typename Number>
class RandomSequence {
    RNG rng;
    std::uniform_real_distribution<Number> dis;
public:
    RandomSequence(RNG&& r, const Number& lower = Number(0), const Number& upper = Number(1)) : 
        rng(std::forward<RNG>(r)), dis(lower,upper) {}

    class const_iterator {
    private:
        RNG rng;
        std::uniform_real_distribution<Number> dis;
        Number n;
        const_iterator(const RNG& r, const std::uniform_real_distribution<Number>& d) :
            rng(r), dis(d), n(dis(rng)) {} 
        friend class RandomSequence<RNG,Number>;
    public:
        const Number& operator*() const { return n; }
        const_iterator& operator++()    {  n = dis(rng); return (*this); }
		bool operator==(const const_iterator& that) const { return false; } //Infinite list
		bool operator!=(const const_iterator& that) const { return true;  } //Infinite list
    }; 

    const_iterator begin() const { return const_iterator(rng,dis); }
    const_iterator end() const { return const_iterator(rng,dis); }
};

template<typename Number,typename RNG = std::mt19937>
auto random_sequence(const Number& lower, const Number& upper, std::size_t seed = std::random_device()()) {
    return RandomSequence<RNG,Number>(RNG(seed),lower,upper);
} 

}
