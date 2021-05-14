#include <Eigen/Dense>

class Material {
    bool emitter_; // True -> emitter. False -> lambertian
    Eigen::Array3f color_;
public:
	Material(bool e, const Eigen::Array3f& c) : emitter_(e), color_(c) { }
	bool is_emitter() const    { return emitter_; }
	bool is_lambertian() const { return !emitter_; }
	const Eigen::Array3f& color() const { return color_; }
};
Material emitter(const Eigen::Array3f& color) {
    return Material(true, color);
}
Material lambertian(const Eigen::Array3f& color) {
    return Material(false, color);
}
#define MATERIAL Material

#include <mj2/tracer/tracer.h>
#include "light-sources.h"

class RenderPathTracing {
    tracer::Scene scene;
    tracer::Pinhole camera;
public:
	template<typename Float>
    Spectrum operator()(const std::array<Float,4>& sample) const {
        auto hit = scene.trace(camera.ray(2.0*sample[0]-1.0f,2.0*sample[1]-1.0f));
        if ( (!hit) || (!hit->material()) ) return Spectrum::Constant(0.0f); 
        else if (hit->material()->is_emitter()) return hit->material()->color();
        else {
            float theta = std::acos(std::min(1.0f,std::sqrt(std::max(0.0f,sample[2]))));
            float phi = 2.0f*M_PI*sample[3];
            tracer::Ray ray(hit->point(),hit->local_to_global()*Eigen::Vector3f(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta)),1.e-3f);
            auto hit2 = scene.trace(ray);
            if ((hit2) && (hit2->material()) && (hit2->material()->is_emitter()))
                return (hit->material()->color()*hit2->material()->color()).eval();
            else
                return Spectrum::Constant(0.0f);
        }
    }

    RenderPathTracing(const tracer::Scene& scene, const tracer::Pinhole& camera) :
        scene(scene), camera(camera) { }
};

template<typename Geometry, typename Hit>
Spectrum direct_light_surface(const Geometry& geometry, const PointLight& light, const Hit& hit) {
	const float eps = 1.e-5;
	Eigen::Vector3f wi = light.position()-hit->point();
	float distance = std::max(wi.norm(),eps);
	wi/=distance;
    return ((!hit)||(!hit->material()))?
		Spectrum::Constant(0.0f).eval():
			(hit->material()->is_emitter()?
				hit->material()->color():
				( ((hit->normal().dot(wi)<=0.0f)||(!geometry.trace_shadow(tracer::Ray(hit->point(),wi,1.e-3f,distance))))?Spectrum::Constant(0.0f).eval():
					hit->material()->color()*light.power(-wi)*hit->normal().dot(wi)/(M_PI*distance*distance)));
}