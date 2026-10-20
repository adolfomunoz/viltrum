#pragma once
#include <mutex>
#include <array>
#include <vector>

namespace viltrum {

template<typename Bins, std::size_t DIMBINS>
class MutexedBins {
    Bins& bins;
    std::array<std::size_t,DIMBINS> bin_resolution;
    std::vector<std::mutex> mutexes;
    using value_type = std::decay_t<decltype(bins(std::array<std::size_t,DIMBINS>{}))>;
    constexpr std::size_t hash_of(const std::array<std::size_t,DIMBINS>& index) {
        std::size_t sum{0}, prod{1};
        for (std::size_t i = 0; i<DIMBINS; ++i) { sum+=prod*index[i]; prod*=bin_resolution[i]; }
        return sum % mutexes.size();
    }
public:
    MutexedBins(Bins& b, const std::array<std::size_t,DIMBINS>& br,std::size_t nmutexes) : 
        bins(b),bin_resolution(br),mutexes(nmutexes) {}
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