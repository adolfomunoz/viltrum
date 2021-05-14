#pragma once
#include "surface.h"

class Medium {
    Spectrum absorption_;
    Spectrum scattering_;
public:
    const Spectrum& absorption() const { return absorption_; }
    const Spectrum& scattering() const { return scattering_; }
    auto extinction() const { return absorption()+scattering();    }
    auto albedo() const { return scattering()/extinction();        }
	Medium(const Spectrum& a, const Spectrum& s) : 
		absorption_(a), scattering_(s) {}
};



template<typename Geometry>
auto medium_light(const Geometry& geometry, const PointLight& light, const tracer::Ray& ray, float t, const Medium& medium) {
    Eigen::Vector3f medium_point = ray.at(t);
//	std::cerr<<"Medium light = "<<((-medium.extinction()*(t + (medium_point - light.position()).norm())).exp()*medium.scattering()*incident_light(geometry,light,medium_point)/(4.0f*M_PI)).maxCoeff()<<" - Incident light = "<<incident_light(geometry,light,medium_point).maxCoeff()<<std::endl;
    return ((-medium.extinction()*(t + (medium_point - light.position()).norm())).exp()*medium.scattering()*incident_light(geometry,light,medium_point)/(4.0f*M_PI)).eval();
}

template<typename Geometry, typename Hit>
Spectrum surface_medium_light(const Geometry& geometry, const PointLight& light, const tracer::Ray& ray, const Hit& hit, const Medium& medium) {
    auto s = direct_light_surface(geometry,light,hit);
    return (s.maxCoeff()<=0.0f)?Spectrum::Constant(0.0f).eval():((-medium.extinction()*((ray.origin() - hit->point()).norm() + (hit->point() - light.position()).norm())).exp()*s);
}

class RenderMediumSingleScattering {
    tracer::Scene scene;
    tracer::Pinhole camera;
    Medium medium;
    PointLight light;
    float max_t;
public:
    Spectrum operator()(const std::array<float,3>& sample) const {
        auto ray = camera.ray(2.0*sample[0]-1.0f,2.0*sample[1]-1.0f);
        auto hit = scene.trace(ray);
        float t = hit?std::min(max_t,hit->distance()):max_t;
//		std::cerr<<"Render function = "<<(surface_medium_light(scene,light,ray,hit,medium) + t*medium_light(scene,light,ray,t*sample[2],medium)).maxCoeff()<<" - t="<<t<<" - Surface = "<<surface_medium_light(scene,light,ray,hit,medium).maxCoeff()<<std::endl;
        return (surface_medium_light(scene,light,ray,hit,medium) + t*medium_light(scene,light,ray,t*sample[2],medium)).eval();
    }

    RenderMediumSingleScattering(const tracer::Scene& scene, const tracer::Pinhole& camera, const Medium& medium, const PointLight& light, float max_t = 10.0f) :
        scene(scene), camera(camera), medium(medium), light(light), max_t(max_t) { }
};

class RenderMediumSingleScatteringDistance {
    tracer::Scene scene;
    tracer::Pinhole camera;
    Medium medium;
    PointLight light;
    float max_t;
public:
	//First value : s, second value factor 
	std::tuple<float,float> sample_distance(float t, float sample) const {
		float sigma_ref = medium.extinction().maxCoeff();
		float e = std::exp(-t*sigma_ref);
		float d = sample*(e-1.0f)+1.0f;
//		if (d<=1.e-5) std::cerr<<sample<<" - "<<e<<" - "<<d<<std::endl;
		return std::tuple<float,float>(-std::log(d)/(t*sigma_ref),-(e-1)/(d*sigma_ref));
	}
	
    Spectrum operator()(const std::array<float,3>& sample) const {
        auto ray = camera.ray(2.0*sample[0]-1.0f,2.0*sample[1]-1.0f);
        auto hit = scene.trace(ray);
        float t = hit?std::min(max_t,hit->distance()):max_t;
		
		auto [s, factor] = this->sample_distance(t, sample[2]);
		
        return (surface_medium_light(scene,light,ray,hit,medium) + 
			factor*medium_light(scene,light,ray,t*s,medium)).eval();
    }

    RenderMediumSingleScatteringDistance(const tracer::Scene& scene, const tracer::Pinhole& camera, const Medium& medium, const PointLight& light, float max_t = 10.0f) :
        scene(scene), camera(camera), medium(medium), light(light), max_t(max_t) { }
};

class RenderMediumSingleScatteringEquiangular {
    tracer::Scene scene;
    tracer::Pinhole camera;
    Medium medium;
    PointLight light;
    float max_t;
public:
	template<typename Hit>
	Spectrum operator()(const tracer::Ray& ray, const Hit& hit, float sample) const {
		const float eps = 1.e-5;
        float t = hit?std::min(max_t,hit->distance()):max_t;
		
		float t_light = (light.position() - ray.origin()).dot(ray.direction());
		//vv Avoid distance 0 so avoid division by 0
		float distance = std::max((light.position() - ray.at(t_light)).norm(),eps);
		
		float r_min = std::atan(-t_light/distance);
		float r = std::atan((t-t_light)/distance) - r_min;
		float s = (distance*std::tan(sample*r + r_min) + t_light)/t;

		Eigen::Vector3f medium_point = ray.at(t*s);
		Eigen::Vector3f d = light.position() - medium_point; 
		float norm = d.norm();
		d/=norm;

		return (surface_medium_light(scene,light,ray,hit,medium) +
			(scene.trace_shadow(tracer::Ray(medium_point, d, 1.e-3f, norm))?
				Spectrum::Constant(0.0f).eval():
				(r/(4.0f*M_PI*distance))*light.power(-d)*(-medium.extinction()*(t + norm)).exp()*medium.scattering())).eval();		
	}
	
    Spectrum operator()(const std::array<float,3>& sample) const {
		auto ray = camera.ray(2.0*sample[0]-1.0f,2.0*sample[1]-1.0f);
		auto hit = scene.trace(ray);
		return (*this)(ray,hit,sample[2]);
    }

    RenderMediumSingleScatteringEquiangular(const tracer::Scene& scene, const tracer::Pinhole& camera, const Medium& medium, const PointLight& light, float max_t = 10.0f) :
        scene(scene), camera(camera), medium(medium), light(light), max_t(max_t) { }
};