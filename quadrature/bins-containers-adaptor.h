#include <vector>
#include "integrate-bins.h"

namespace viltrum {
    
template <typename T>
class adaptor_vector {
    std::vector<T>& bins;
public:
    adaptor_vector(std::vector<T>& v) : bins(v) {}
    std::array<std::size_t,1> resolution() const { return std::array<std::size_t,1>{bins.size()}; }
    T& operator()(const std::array<std::size_t,1>& i) { return bins[i[0]]; }
    const T& operator()(const std::array<std::size_t,1>& i) const { return bins[i[0]]; }
}; 

template<typename T>
adaptor_vector<T> adaptor(std::vector<T>& v) {
    return adaptor_vector<T>(v);
}   

template <typename T>
class adaptor_vector_2d {
    std::vector<std::vector<T>>& bins;
public:
    adaptor_vector_2d(std::vector<std::vector<T>>& v) : bins(v) {
        //Adjust size to max row size so it is regular
        if (bins.size()<1) bins.resize(1);
        std::size_t s = 1;
        for (auto& row : bins) if (row.size()>s) s=row.size();
        for (auto& row : bins) if (row.size()!=s) row.resize(s);        
    }
    std::array<std::size_t,2> resolution() const { return std::array<std::size_t,2>{bins[0].size(),bins.size()}; }
    T& operator()(const std::array<std::size_t,2>& i) { return bins[i[1]][i[0]]; }
    const T& operator()(const std::array<std::size_t,2>& i) const { return bins[i[1]][i[0]]; }
}; 

template<typename T>
adaptor_vector_2d<T> adaptor(std::vector<std::vector<T>>& v) {
    return adaptor_vector_2d<T>(v);
}    
 
    
template<typename IntegratorBins, typename Container, typename F, typename Float, std::size_t DIM>
void integrate_bins(const IntegratorBins& integrator_bins, Container& bins, const F& function, const Range<Float,DIM>& range) {
    auto a = adaptor(bins);
    integrate_bins(integrator_bins, a, a.resolution(), function, range);
}

    
}