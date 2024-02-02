#pragma once

#include <vector>
#include <array>

namespace viltrum {

template<typename T, std::size_t DIMBINS>
class tensor {
    std::vector<T> data;
    std::array<std::size_t, DIMBINS> res;
public:
    constexpr const std::array<std::size_t, DIMBINS>& resolution() const { return res; }
    constexpr std::size_t resolution(std::size_t i) const { return resolution()[i]; }

private:
    std::size_t position(const std::array<std::size_t,DIMBINS>& p) const {
        std::size_t pos = 0, prod = 1;
        for (std::size_t d = 0;d<DIMBINS; ++d) {
            pos += p[d]*prod; prod*=resolution(d);
        }
        return pos;
    }
public:
    tensor(const std::array<std::size_t, DIMBINS>& r, const T& t = T()) : res(r) {
        std::size_t elements(1);
        for (auto r : res) elements*=r;
        data.resize(elements, t);
    }
    
    tensor(const std::array<std::size_t, DIMBINS>& r, const std::vector<T>& e) : res(r),data(e) {
        std::size_t elements(1);
        for (auto r : res) elements*=r;
        data.resize(elements);
    }
    
    const std::vector<T>& raw_data() const { return data; } 

    T& operator[](const std::array<std::size_t,DIMBINS>& p) {
        return data[position(p)];
    }

    T& operator()(const std::array<std::size_t,DIMBINS>& p) {
        return (*this)[p];
    }

    const T& operator[](const std::array<std::size_t,DIMBINS>& p) const {
        return data[position(p)];
    }

    const T& operator()(const std::array<std::size_t,DIMBINS>& p) const {
        return (*this)[p];
    }


    std::size_t size() const { return raw_data().size(); }
};

}

