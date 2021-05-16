#pragma once

#include <random>

namespace viltrum {

template<typename RNG>
class VectorSamplerUniform {
	mutable RNG rng;
		
	class Sampler {
		mutable RNG rng;
		std::size_t size;
	public:
		//Returns position and probability
		std::tuple<std::size_t,double> sample() {
			std::uniform_int_distribution<std::size_t> choose(0,size-1);
			return std::make_tuple(choose(rng),1.0/double(size));
		}
		
		Sampler(const RNG& r, std::size_t s) : rng(r), size(s) {}
		Sampler() {} //This makes things easier although I don't like it
	};

public:
	template<typename T>
	Sampler operator()(const std::vector<T>& v) const {
		std::uniform_int_distribution<std::size_t> choose(0,10000000);
		return Sampler(RNG(choose(rng)),v.size());
	}

    VectorSamplerUniform(RNG&& r) : rng(std::forward<RNG>(r)) { }
};

template<typename RNG>
auto vector_sampler_uniform(RNG&& rng) {
    return VectorSamplerUniform<RNG>(std::forward<RNG>(rng));
}

auto vector_sampler_uniform(std::size_t seed = std::random_device()()) {
    return vector_sampler_uniform(std::mt19937_64(seed));
}

}



