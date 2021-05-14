#pragma once
#include <Eigen/Dense>
struct Material {
    bool emitter; // True -> emitter. False -> lambertian
    Eigen::Array3f color;
};
#define MATERIAL Material
#include <mj2/tracer/tracer.h>
#include <functional>
#include "../quadrature/monte-carlo.h"
#include <PerlinNoise/PerlinNoise.hpp>
#include <cimg-all.h>

Material emitter(const Eigen::Array3f& color) {
    return Material{true, color};
}
Material lambertian(const Eigen::Array3f& color) {
    return Material{false, color};
}


class RenderAreaLighting {
    tracer::Scene scene;
    tracer::Pinhole camera;
public:
    Eigen::Array3f operator()(const std::array<float,4>& sample) const {
        auto hit = scene.trace(camera.ray(2.0*sample[0]-1.0f,2.0*sample[1]-1.0f));
        if ( (!hit) || (!hit->material()) ) return Eigen::Array3f(1.e-5,1.e-5,1.e-5); 
        else if (hit->material()->emitter) return hit->material()->color;
        else {
            float theta = std::acos(std::min(1.0f,std::sqrt(std::max(0.0f,sample[2]))));
            float phi = 2.0f*M_PI*sample[3];
            tracer::Ray ray(hit->point(),hit->local_to_global()*Eigen::Vector3f(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta)),1.e-3f);
            auto hit2 = scene.trace(ray);
            if ((hit2) && (hit2->material()) && (hit2->material()->emitter))
                return hit->material()->color*hit2->material()->color;
            else
                return Eigen::Array3f(1.e-5,1.e-5,1.e-5);
        }
    }

    RenderAreaLighting(const tracer::Scene& scene, const tracer::Pinhole& camera) :
        scene(scene), camera(camera) { }
};

tracer::Scene& environment_constant(tracer::Scene& scene) {
    scene.add(tracer::Sphere(Eigen::Vector3f(0,0,0),0.5).set_material(lambertian(Eigen::Array3f(0.8,0.8,0.8))));
    scene.add(tracer::Sphere(Eigen::Vector3f(0,0,0),5).set_material(emitter(Eigen::Array3f(1,1,1))));
    scene.add(tracer::Plane(Eigen::Vector3f(0,1,0),Eigen::Vector3f(0,-0.5,0)).set_material(lambertian(Eigen::Array3f(0.8,0.8,0.8))));
	return scene;
}

std::tuple<std::function<Eigen::Array3f(const std::array<float,4>&)>,std::function<Eigen::Array3f(const std::array<float,4>&,const std::array<float,4>&)>> render_function4d(int argc, char **argv) {
	std::function<Eigen::Array3f(const std::array<float,4>&)> func = [] (const std::array<float,4>& x) { return Eigen::Array3f(1.0f,1.0f,1.0f); };

	std::function<Eigen::Array3f(const std::array<float,4>&,const std::array<float,4>&)> groundtruth;

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
			tracer::Pinhole camera(Eigen::Vector3f( 0, 0, -3), Eigen::Vector3f( 0, 0, 2), Eigen::Vector3f( 0, 1, 0), Eigen::Vector3f(w/h, 0, 0));
			tracer::Scene scene;
			if (std::string(argv[i])=="area-constant") {
				environment_constant(scene);
			}
			func = RenderAreaLighting(scene,camera);
		}
	}
    if (!has_groundtruth) groundtruth = [func,ground_truth_samples,seed] (const std::array<float,4>& a, const std::array<float,4>& b) {
        return integrator_monte_carlo_uniform(ground_truth_samples,seed).integrate(func,range(a,b));
	};

	return std::make_tuple(func,groundtruth);
}

