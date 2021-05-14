#pragma once
#include "medium.h"

class RenderMediumTwoBounces {
	RenderMediumSingleScatteringEquiangular equiangular;
	RenderMediumSingleScatteringDistance distance;
	tracer::Scene scene;
    tracer::Pinhole camera;
	Medium medium;
	float max_t;
public:
    Spectrum operator()(const std::array<float,6>& sample) const {
        auto ray = camera.ray(2.0*sample[0]-1.0f,2.0*sample[1]-1.0f);
        auto hit = scene.trace(ray);
        float t = hit?std::min(max_t,hit->distance()):max_t;
		
		float s, factor;
		std::tie(s, factor) = distance.sample_distance(t, sample[2]);
		
//		std::cerr<<"Medium hit at "<<s*t<<" out of "<<t<<std::endl;
		float phi = 2.0f*M_PI*sample[3];
		float cos_theta = std::max(-1.0f,std::min(1.0f,2.0f*sample[4]-1.0f));
		float sin_theta = std::sqrt(1.0f - cos_theta*cos_theta);
//		std::cerr<<"Secondary direction : "<<sin_theta*std::cos(phi)<<", "<<sin_theta*std::sin(phi)<<", "<<cos_theta<<std::endl;
		tracer::Ray secondary_ray_medium(ray.at(s*t),Eigen::Vector3f(sin_theta*std::cos(phi),sin_theta*std::sin(phi), cos_theta));


		Spectrum surface_indirect(0);
		if ((hit) && (hit->material()) && (hit->material()->is_lambertian())) {
			float theta = std::acos(std::min(1.0f,std::sqrt(std::max(0.0f,sample[3]))));
            float phi = 2.0f*M_PI*sample[4];
            tracer::Ray secondary_ray_surface(hit->point(),hit->local_to_global()*Eigen::Vector3f(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta)),1.e-3f);
            auto hit2 = scene.trace(secondary_ray_surface);
			surface_indirect = ((-medium.extinction()*hit->distance()).exp()*hit->material()->color()*equiangular(secondary_ray_surface,hit2,sample[5])).eval();
		}
		
/*		if (((equiangular(ray,hit,sample[2]) + (factor*(-medium.extinction()*(s*t)).exp()*medium.scattering())*equiangular(secondary_ray_medium,scene.trace(secondary_ray_medium),sample[5]) + surface_indirect).maxCoeff())>1.e4) {
			std::cerr<<equiangular(ray,hit,sample[2]).maxCoeff()<<" - ";
			std::cerr<<((factor*(-medium.extinction()*(s*t)).exp()*medium.scattering())*equiangular(secondary_ray_medium,scene.trace(secondary_ray_medium),sample[5]) + surface_indirect).maxCoeff()<<" - ";
			std::cerr<<surface_indirect.maxCoeff()<<std::endl;
		}
		*/
		
		return (equiangular(ray,hit,sample[2]) + (factor*(-medium.extinction()*(s*t)).exp()*medium.scattering())*equiangular(secondary_ray_medium,scene.trace(secondary_ray_medium),sample[5]) + surface_indirect).eval();
//		return ((factor*(-medium.extinction()*(s*t)).exp()*medium.scattering())*equiangular(secondary_ray_medium,scene.trace(secondary_ray_medium),sample[5])).eval();
    }
	
	 RenderMediumTwoBounces(const tracer::Scene& scene, const tracer::Pinhole& camera, const Medium& medium, const PointLight& light, float max_t = 10.0f) : equiangular(scene,camera,medium,light,max_t), distance(scene,camera,medium,light,max_t), scene(scene),camera(camera),medium(medium),max_t(max_t) {}
};
/*
class RenderMediumSecoundBounce {
	RenderMediumSingleScatteringEquiangular equiangular;
	RenderMediumSingleScatteringDistance distance;
	tracer::Scene scene;
    tracer::Pinhole camera;
	Medium medium;
	float max_t;
public:
    Spectrum operator()(const std::array<float,6>& sample) const {
        auto ray = camera.ray(2.0*sample[0]-1.0f,2.0*sample[1]-1.0f);
        auto hit = scene.trace(ray);
        float t = hit?std::min(max_t,hit->distance()):max_t;
		
		auto [s, factor] = distance.sample_distance(t, sample[2]);
		
		float phi = 2.0f*M_PI*sample[3];
		float cos_theta = 2.0f*sample[4]-1.0f;
		float theta = std::acos(cos_theta);
		tracer::Ray secondary_ray_medium(ray.at(s*t),Eigen::Vector3f(std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi), cos_theta));
		
		return ((-medium.extinction()*(s*t)).exp()*medium.scattering()*factor*equiangular(secondary_ray_medium,scene.trace(secondary_ray_medium),sample[5])).eval();
    }
	
	 RenderMediumSecondBounce(const tracer::Scene& scene, const tracer::Pinhole& camera, const Medium& medium, const PointLight& light, float max_t = 10.0f) : equiangular(scene,camera,medium,light,max_t), distance(scene,camera,medium,light,max_t), scene(scene),camera(camera),medium(medium),max_t(max_t) {}
};

class RenderMediumSurface {
	RenderMediumSingleScatteringEquiangular equiangular;
	tracer::Scene scene;
    tracer::Pinhole camera;
	Medium medium;
	float max_t;
public:
    Spectrum operator()(const std::array<float,6>& sample) const {
        auto ray = camera.ray(2.0*sample[0]-1.0f,2.0*sample[1]-1.0f);
        auto hit = scene.trace(ray);
 
		if ((hit) && (hit->material()) && (hit->material()->is_lambertian())) {
			float theta = std::acos(std::min(1.0f,std::sqrt(std::max(0.0f,sample[3]))));
            float phi = 2.0f*M_PI*sample[4];
            tracer::Ray secondary_ray_surface(hit->point(),hit->local_to_global()*Eigen::Vector3f(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta)),1.e-3f);
            auto hit2 = scene.trace(secondary_ray_surface);
			return ((-medium.extinction()*hit->distance()).exp()*hit->material()->color()*equiangular(secondary_ray_surface,scene.trace(secondary_ray_surface),sample[5])).eval();
		} else return Spectrum::Constant(0).eval();
    }
	
	 RenderMediumSurface(const tracer::Scene& scene, const tracer::Pinhole& camera, const Medium& medium, const PointLight& light, float max_t = 10.0f) : equiangular(scene,camera,medium,light,max_t), scene(scene),camera(camera),medium(medium),max_t(max_t) {}
};
*/