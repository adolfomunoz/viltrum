#pragma once
#include "nested.h"
#include "error-heuristic.h"


namespace viltrum {

template<typename Rule, typename ErrorHeuristic, typename = std::enable_if_t<is_nested<Rule>::value>>
class RegionsGeneratorAdaptiveHeap {
    Rule rule;
    ErrorHeuristic error_heuristic;
    std::size_t subdivisions;

public:
    RegionsGeneratorAdaptiveHeap(const Rule& r, const ErrorHeuristic& er, std::size_t subdivisions) 
        : rule(r), error_heuristic(er), subdivisions(subdivisions) {}

    template<std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
    auto generate(const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {
            auto heap_ordering = [] (const auto& a, const auto& b) {
                return std::get<0>(a.extra()) < std::get<0>(b.extra());
            };


            logger.log_progress(std::size_t(0),subdivisions);
            auto r = region(f,rule,range.min(),range.max());
            auto errdim = error_heuristic(r);
            std::vector<ExtendedRegion<decltype(r),decltype(errdim)> > heap;
            heap.reserve(subdivisions+1);
            heap.emplace_back(r,errdim);
            for (std::size_t i = 0; i<subdivisions;++i) {
    	        auto r = heap.front();
	            auto subregions = r.split(f,std::get<1>(r.extra()));
	            std::pop_heap(heap.begin(),heap.end(),heap_ordering); heap.pop_back();
	            for (auto sr : subregions) {
		            auto errdim = error_heuristic(sr);
		            heap.emplace_back(sr,errdim); 
		            std::push_heap(heap.begin(), heap.end(), heap_ordering);
                }
                logger.log_progress(i,subdivisions);
            }
            logger.log_progress(subdivisions,subdivisions);
            return heap;
	    }
};

template<typename Rule, typename ErrorHeuristic>
auto regions_generator_adaptive_heap(const Rule& r, const ErrorHeuristic& er, std::size_t subdivisions) {
    return RegionsGeneratorAdaptiveHeap<Rule,ErrorHeuristic>(r,er,subdivisions);
}


}