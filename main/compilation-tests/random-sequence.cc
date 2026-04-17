#include "../../src/monte-carlo/random-sequence-ref.h"
#include "../../src/monte-carlo/random-sequence-rng.h"
#include "../../src/monte-carlo/random-sequence-seed.h"
#include "../../src/rng/xoroshiro.h"
#include <iostream>
#include <chrono>



using namespace viltrum;


template<typename Seq>
void test_sequence(const char* name, const Seq& seq, unsigned long n = 10000000) {
    unsigned long i = 0;
    float sum = 0.0f;

    auto t0 = std::chrono::high_resolution_clock::now();
    for (float x : seq) {
        sum += x;
        if (++i >= n) break;
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::cout << name<<"  - "<<sum << "  (elapsed: " << elapsed << " s)" << std::endl;
}

int main() {
    RangeInfinite<float> range(std::vector{-2.0f},std::vector{2.0f});
    test_sequence("Seed - mt19937",random_sequence_seed(range,0));
    test_sequence("RNG  - mt19937",random_sequence_rng(range,0));
    std::mt19937 rng(0);
    test_sequence("Ref  - mt19937",random_sequence_ref(range,rng));
    
    test_sequence("Seed - xoros. ",random_sequence_seed<float,Xoroshiro128Plus>(range,0));
    test_sequence("RNG  - xoros. ",random_sequence_rng<float,Xoroshiro128Plus>(range,0));
    Xoroshiro128Plus rng2(0);
    test_sequence("Ref  - xoros. ",random_sequence_ref(range,rng2));
    std::cout<<std::endl;
}  