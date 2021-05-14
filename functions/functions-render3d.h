#pragma once
#include <Eigen/Dense>
struct Material {
    Eigen::Array3f color; //Lambertian diffuse surfaces (all of them)
};
#define MATERIAL Material
#include <mj2/tracer/tracer.h>
#include <functional>
#include "../quadrature/monte-carlo.h"
#include <PerlinNoise/PerlinNoise.hpp>
#include <cimg-all.h>

Material lambertian(const Eigen::Array3f& color) {
    return Material{color};
}

struct Medium {
    Eigen::Array3f absorption_;
    Eigen::Array3f scattering_;
    const Eigen::Array3f& absorption() const { return absorption_; }
    const Eigen::Array3f& scattering() const { return scattering_; }
    auto extinction() const { return absorption()+scattering();    }
    auto albedo() const { return scattering()/extinction();        }
};

struct PointLight {
    Eigen::Array3f power_;
    Eigen::Vector3f position_;
    const Eigen::Array3f& power() const { return power_; }
    const Eigen::Vector3f& position() const { return position_; }
    const Eigen::Array3f& power_towards(const Eigen::Vector3f& d) const { return power(); }
};

struct ConeLight {
    Eigen::Array3f power_;
    Eigen::Vector3f position_;
    Eigen::Vector3f direction_;
    float angle_;
    const Eigen::Array3f& power() const { return power_; }
    const Eigen::Vector3f& position() const { return position_; }
    Eigen::Array3f power_towards(const Eigen::Vector3f& d) const {
        return ((d.dot(direction_)/d.norm())>=std::cos(angle_))?power():Eigen::Array3f(0.0f,0.0f,0.0f);
    }
};


template<typename Light>
class RenderMediumSingleScattering {
    tracer::Scene scene;
    tracer::Pinhole camera;
    Medium medium;
    Light light;
    float max_t;
public:
    Eigen::Array3f operator()(const std::array<float,3>& sample) const {
        float t; Eigen::Array3f back;
        auto ray = camera.ray(2.0*sample[0]-1.0f,2.0*sample[1]-1.0f);
        auto hit = scene.trace(ray);
        if (!hit) {
            t = max_t; back = Eigen::Array3f(0,0,0);
        } else { 
            t = std::min(max_t,hit->distance());
            if (!hit->material()) back = Eigen::Array3f(0,0,0); 
            else if ( (hit->normal().dot(light.position() - hit->point())<=0.0f) || (scene.trace_shadow(tracer::Ray(hit->point(),light.position()-hit->point(),1.e-3f,(light.position()-hit->point()).norm()))) ) back = Eigen::Array3f(0,0,0);
            else back = hit->material()->color*light.power_towards(hit->point() - light.position())*hit->normal().dot(light.position() - hit->point())/(M_PI*std::pow((hit->point() - light.position()).norm(),3));
        }

        Eigen::Vector3f medium_point = ray.at(t*sample[2]);
        Eigen::Array3f sol = (-medium.extinction()*t).exp()*back;
        if (!scene.trace_shadow(tracer::Ray(medium_point,light.position()-medium_point,1.e-3f,(light.position()-medium_point).norm()))) {
            sol += t*(-medium.extinction()*(t*sample[2] + (medium_point - light.position()).norm())).exp()*medium.scattering()*light.power_towards(medium_point - light.position())/(4.0f*M_PI);
        }
        return sol;
    }

    RenderMediumSingleScattering(const tracer::Scene& scene, const tracer::Pinhole& camera, const Medium& medium, const Light& light, float max_t = 10.0f) :
        scene(scene), camera(camera), medium(medium), light(light), max_t(max_t) { }
};


std::tuple<std::function<Eigen::Array3f(const std::array<float,3>&)>,std::function<Eigen::Array3f(const std::array<float,3>&,const std::array<float,3>&)>> render_function3d(int argc, char **argv) {
	std::function<Eigen::Array3f(const std::array<float,3>&)> func = [] (const std::array<float,3>& x) { return Eigen::Array3f(1.0f,1.0f,1.0f); };

	std::function<Eigen::Array3f(const std::array<float,3>&,const std::array<float,3>&)> groundtruth;

    unsigned long ground_truth_samples = 10000000;
    std::size_t seed = std::random_device()();
    bool has_groundtruth = false;
	int w = 512; int h=512;
    for (int i = 0; i < argc-1; ++i) {
		if (std::string(argv[i])=="-width") w = atoi(argv[++i]);
		else if (std::string(argv[i])=="-height") h = atoi(argv[++i]);
        if (std::string(argv[i])=="-ground-truth-samples") ground_truth_samples = atol(argv[++i]);
        else if (std::string(argv[i])=="-ground-truth-seed") seed = atol(argv[++i]);
    }

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-scene") { ++i;
			if (std::string(argv[i])=="single-scattering") {
                ++i;
			    tracer::Pinhole camera(Eigen::Vector3f( 0, 0, -3), Eigen::Vector3f( 0, 0, 2), Eigen::Vector3f( 0, 1, 0), Eigen::Vector3f(w/h, 0, 0));
//                ConeLight light{Eigen::Array3f(1,1,1),Eigen::Vector3f(2,0,0),Eigen::Vector3f(-1,0,0),M_PI/8.0f};;
                PointLight light{Eigen::Array3f(1,1,1),Eigen::Vector3f(2,0,0)};;
			    tracer::Scene scene;
                if ((i<argc) && (std::string(argv[i])=="occluded")) { ++i; scene.add(tracer::Sphere(Eigen::Vector3f(0.5,0,0),0.25).set_material(lambertian(Eigen::Array3f(0.2,0.8,0.2)))); }
                scene.add(tracer::Plane(Eigen::Vector3f(1,0,0),Eigen::Vector3f(-1,0,0)).set_material(lambertian(Eigen::Array3f(0.8,0.2,0.2))));
			    func = RenderMediumSingleScattering(scene,camera,Medium{Eigen::Array3f(0.05f,0.05f,0.05f),Eigen::Array3f(0.25f,0.25f,0.25f)},light,6.0f);
			}
		}
	}
    if (!has_groundtruth) groundtruth = [func,ground_truth_samples,seed] (const std::array<float,3>& a, const std::array<float,3>& b) {
        return integrator_monte_carlo_uniform(ground_truth_samples,seed).integrate(func,range(a,b));
	};

	return std::make_tuple(func,groundtruth);
}

