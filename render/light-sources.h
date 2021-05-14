#pragma once
#include <Eigen/Dense>
#include <memory>
#include <cmath>
#include <type_traits>


#ifndef Spectrum
using Spectrum = Eigen::Array3f;
#endif

class SphericalSpectrumBase {
public:
	virtual Spectrum value(const Eigen::Vector3f& d) const = 0;	
};

class SphericalSpectrumConstant : public SphericalSpectrumBase {
	Spectrum s;
public:
	Spectrum value(const Eigen::Vector3f& d) const override { return s; }
	SphericalSpectrumConstant(const Spectrum& s) : s(s) {}
};

class SphericalSpectrumCone : public SphericalSpectrumBase {
	Eigen::Vector3f direction; float angle;
	Spectrum s;
public:
	Spectrum value(const Eigen::Vector3f& d) const override {
		return (direction.dot(d)>=std::cos(angle))?s:Spectrum::Constant(0.0f);
	}
	
	SphericalSpectrumCone(const Eigen::Vector3f& direction, float angle,
		const Spectrum& s) :
			direction(direction), angle(angle), s(s) {}
};
	
class SphericalSpectrum {
	std::shared_ptr<SphericalSpectrumBase> spectrum;
public:
	SphericalSpectrum(const Spectrum& s) :
		spectrum(std::make_shared<SphericalSpectrumConstant>(s)) {}
	template<typename SS>
	SphericalSpectrum(const SS& ss, std::enable_if_t<std::is_base_of_v<SphericalSpectrumBase,SS>>* p = nullptr) :
		spectrum(std::make_shared<std::decay_t<SS>>(ss)) {}
	template<typename SS>
	SphericalSpectrum(SS&& ss, std::enable_if_t<std::is_base_of_v<SphericalSpectrumBase,SS>>* p = nullptr) :
		spectrum(std::make_shared<std::decay_t<SS>>(std::forward<SS>(ss))) {}
		
	Spectrum operator()(const Eigen::Vector3f& d) const { return spectrum->value(d); }
};

class PointLight {
    Eigen::Vector3f position_;
public:
	SphericalSpectrum power;
    const Eigen::Vector3f& position() const { return position_; }
    Spectrum power_at(const Eigen::Vector3f& point) const {
		Eigen::Vector3f direction = point - position();
		float f = 1.0f/direction.squaredNorm();
		direction /= direction.norm();
		return f*power(direction); 
	}
	template<typename... Rest>
	PointLight(const Eigen::Vector3f& position, Rest&&... rest) :
		position_(position), power(std::forward<Rest>(rest)...) {}
};

template<typename Geometry>
auto incident_light(const Geometry& geometry, const PointLight& light, const Eigen::Vector3f& point) {
    Eigen::Vector3f d = light.position() - point;
	float norm = d.norm();
	d/=norm;
    return geometry.trace_shadow(tracer::Ray(point, d, 1.e-3f, norm))?Spectrum::Constant(0.0f).eval():light.power(-d)/(norm*norm);
}
