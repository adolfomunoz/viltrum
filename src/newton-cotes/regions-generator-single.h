#pragma once

namespace viltrum {

template<typename Rule>
class RegionsGeneratorSingle {
    Rule rule;
public:
    RegionsGeneratorSingle(const Rule& r) : rule(r) {}

    template<std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
    auto generate(const std::array<std::size_t,DIMBINS>& bin_resolution,
		const F& f, const Range<Float,DIM>& range, Logger& logger) const {
            logger.log_progress(0,1);
            auto r = region(f,rule, range.min(), range.max());
            logger.log_progress(1,1);
            return std::vector<decltype(r)>(1,r);
        }
};

template<typename Rule>
auto regions_generator_single(const Rule& r) {
    return RegionsGeneratorSingle<Rule>(r);
}


}