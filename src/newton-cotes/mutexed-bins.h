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
};

}