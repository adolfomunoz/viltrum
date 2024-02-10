#pragma once
#include <mutex>
#include <array>
#include <vector>
#include "../tensor.h"

namespace viltrum {
template<typename T, std::size_t DIMBINS>
class MutexedTensorVector {
    tensor<std::vector<T>,DIMBINS> data;
    std::vector<std::mutex> mutexes;
    constexpr std::size_t hash_of(const std::array<std::size_t,DIMBINS>& index) {
        std::size_t sum{0}, prod{1};
        for (std::size_t i = 0; i<DIMBINS; ++i) { sum+=prod*index[i]; prod*=data.resolution(i); }
        return sum % mutexes.size();
    }
public:
    MutexedTensorVector(const std::array<std::size_t, DIMBINS>& r,std::size_t nmutexes) : 
        data(r), mutexes(nmutexes) {}

    void push_at(const std::array<std::size_t,DIMBINS>& i,T&& v) {
        std::scoped_lock lock(mutexes[hash_of(i)]);
        data[i].push_back(std::forward<T>(v));
    }
    //Read only access
    const std::vector<T>& operator[](const std::array<std::size_t,DIMBINS>& p) const {
        return data[p];
    }
};

}