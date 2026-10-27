#pragma once
#include <array>
#include <atomic>
#include <cstddef>
#include <mutex>
#include <thread>
#include <type_traits>
#include <vector>

namespace viltrum {

/**
 * This has been partially vibe-coded.
 * 
 * Same interface as MutexedBins, but more efficient:
 *  - Each lock lives in its own cache line, so threads locking different bins never
 *    invalidate each other's lock (no false sharing).
 *  - The lock is a lightweight test-and-test-and-set spinlock (1 byte of state), since the 
 *    critical section is just a single accumulation.
 *  - The number of locks is rounded up to a power of two so the hash is a mask instead of a modulo,
 *    and the strides of the bin index are precomputed.
 */
template<typename Bins, std::size_t DIMBINS>
class MutexedBinsAligned {
    static constexpr std::size_t cache_line = 64;

    struct alignas(cache_line) Lock {
        std::atomic<bool> locked{false};
        void lock() noexcept {
            unsigned spins = 0;
            while (locked.exchange(true,std::memory_order_acquire)) {
                while (locked.load(std::memory_order_relaxed)) {
                    if (++spins > 64) { std::this_thread::yield(); spins = 0; }
                }
            }
        }
        void unlock() noexcept { locked.store(false,std::memory_order_release); }
    };

    Bins& bins;
    std::array<std::size_t,DIMBINS> strides;
    std::size_t mask;
    std::vector<Lock> mutexes;
    using value_type = std::decay_t<decltype(bins(std::array<std::size_t,DIMBINS>{}))>;

    static std::size_t round_up_pow2(std::size_t n) { 
        std::size_t p = 1;
        while (p < n) p <<= 1;
        return p;
    }

    std::size_t hash_of(const std::array<std::size_t,DIMBINS>& index) const noexcept {
        std::size_t sum{0};
        for (std::size_t i = 0; i<DIMBINS; ++i) sum+=strides[i]*index[i];
        return sum & mask;
    }
public:
    MutexedBinsAligned(Bins& b, const std::array<std::size_t,DIMBINS>& br,std::size_t nmutexes) : 
        bins(b),mask(round_up_pow2(nmutexes)-1),mutexes(mask+1) {
        std::size_t prod{1};
        for (std::size_t i = 0; i<DIMBINS; ++i) { strides[i] = prod; prod*=br[i]; }
    }
    void add(const std::array<std::size_t,DIMBINS>& i,const value_type& v) {
        std::scoped_lock lock(mutexes[hash_of(i)]);
        bins(i) += v;
    }
    // This is a template version of the add function that allows adding values of different types to the bins if 
    // the V can be added to bins' value_type. This is useful for cases where the bins are of a certain type (e.g., double)
    // but we add another type (e.g., float) that can be implicitly converted to the bins' value_type.
    template<typename V>
    void add(const std::array<std::size_t,DIMBINS>& i,const V& v) {
        std::scoped_lock lock(mutexes[hash_of(i)]);
        bins(i) += v;
    }
};

}
