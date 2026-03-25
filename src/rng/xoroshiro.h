#include <cstdint>
#include <limits>
#include <random>
#include <iostream>

namespace viltrum {

class Xoroshiro128Plus {
public:
    // 1. C++11 URBG standard requirements
    using result_type = uint64_t;

    static constexpr result_type min() {
        return std::numeric_limits<result_type>::min();
    }

    static constexpr result_type max() {
        return std::numeric_limits<result_type>::max();
    }

    // 2. Constructors
    explicit Xoroshiro128Plus(uint64_t seed_val = 1) {
        seed(seed_val);
    }

    // 3. Seeding function
    void seed(uint64_t seed_val) {
        // The authors of Xoroshiro strictly recommend using SplitMix64 
        // to initialize the state. Xoroshiro's state can never be all zeros, 
        // and SplitMix64 cleanly mixes a 64-bit seed into 128 bits of state.
        s[0] = splitmix64(seed_val);
        s[1] = splitmix64(seed_val);
    }

    // 4. The generator function
    result_type operator()() {
        const uint64_t s0 = s[0];
        uint64_t s1 = s[1];
        const uint64_t result = s0 + s1;

        s1 ^= s0;
        s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
        s[1] = rotl(s1, 37); // c

        return result;
    }

private:
    uint64_t s[2];

    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    // SplitMix64 used exclusively for seeding
    static uint64_t splitmix64(uint64_t& x) {
        uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
};

}