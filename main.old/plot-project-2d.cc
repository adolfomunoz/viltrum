#include "../viltrum.h"
#include "../functions/functions2d.h"
#include "../plot/integration2d.h"
#include <svg-cpp-plot/svg-cpp-plot.h>

#include <iostream>
#include <cmath>
#include <chrono>

using namespace viltrum;

template<typename F>
class FunctionWrapper {
	F f;
	mutable std::shared_ptr<unsigned long> nsamples;
public:
	FunctionWrapper(const F& f) : f(f), 
		nsamples(std::make_shared<unsigned long>(0)) {}
	FunctionWrapper(F&& f) : f(std::forward<F>(f)), 
		nsamples(std::make_shared<unsigned long>(0)) {}
	
	double operator()(const std::array<double,2>& pos) const {
		++(*nsamples);
		return std::apply(f,pos);		
	}
		
	unsigned long samples() const { return *nsamples; }
};

template<typename IntegratorBins, typename F>
svg_cpp_plot::_2d::Group plot(const char* name, const IntegratorBins& integrator, const F& function, const Range<double,2>& range, unsigned int bins, int resolution = 200, const std::vector<double>& gt = std::vector<double>()) {
    auto sol = svg_cpp_plot::_2d::group();
    auto& graph_function = sol.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,-210)))
                .add(svg_cpp_plot::Graph2D({200,200},svg_cpp_plot::BoundingBox(range.min(0),range.min(1),range.max(0),range.max(1))));
    graph_function.area().add(svg_cpp_plot::_2d::function_2d(function,svg_cpp_plot::_2d::color_map_heat(0.0f,1.0f),{range.min(0),range.min(1)},{range.max(0),range.max(1)},{resolution,resolution}));
    graph_function.border().stroke_width(1).stroke(svg_cpp_plot::black);

    FunctionWrapper f(function);
    std::vector<double> result(bins,0);
    auto vb = vector_bins(result);
	auto start = std::chrono::steady_clock::now();
    integrator.integrate(vb, vector_resolution(result), f, range);
    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);    
	auto& projected = sol.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,10))).add(svg_cpp_plot::Graph2D({200,100},svg_cpp_plot::BoundingBox(0,0,1,1)));
    projected.border().stroke_width(1).stroke(svg_cpp_plot::black);
	projected.area().add(svg_cpp_plot::_2d::bar_plot(result)).stroke_width(0.5).stroke(svg_cpp_plot::black).fill(svg_cpp_plot::yellow);
	if (gt.size()>0)
		projected.area().add(svg_cpp_plot::_2d::bar_plot(gt)).stroke_width(0.5).stroke(svg_cpp_plot::red).fill(svg_cpp_plot::none);

	std::stringstream label_samples;
	label_samples<<std::fixed<<f.samples()<<" samples";
	std::stringstream label_time;
	label_time<<std::fixed<<std::setprecision(3)<<elapsed.count()<<" seconds";

	sol.add(svg_cpp_plot::_2d::text({100,120},label_samples.str())).font_size(10).text_anchor(svg_cpp_plot::text_anchor_middle);
	sol.add(svg_cpp_plot::_2d::text({100,130},label_time.str())).font_size(10).text_anchor(svg_cpp_plot::text_anchor_middle);

    sol.add(svg_cpp_plot::_2d::text({100,-218},std::string(name))).font_size(18).text_anchor(svg_cpp_plot::text_anchor_middle);
	
	graph_function.area().add(plot_samples_2d(integrator,f,range,bins)).stroke_width(2).stroke(svg_cpp_plot::rgb(0.7,0.7,0.7));

    return sol;
}

template<typename Nested, typename Error, typename F>
svg_cpp_plot::_2d::Group plot_adaptive(Nested&& nested, Error&& error, const F& function, const Range<double,2>& range, unsigned long iterations) {
	FunctionWrapper f(function);
    auto sol = svg_cpp_plot::_2d::group();
    auto& graph_function = sol.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,-210)))
                .add(svg_cpp_plot::Graph2D({200,200},svg_cpp_plot::BoundingBox(range.min(0),range.min(1),range.max(0),range.max(1))));
	graph_function.area().add(plot_adaptive_boundaries_2d(nested,error,f,range,iterations)).stroke_width(1).stroke(svg_cpp_plot::rgb(1,1,1));
	graph_function.area().add(plot_adaptive_samples_2d(nested,error,f,range,iterations)).stroke_width(2).stroke(svg_cpp_plot::rgb(1,1,1));
    return sol;
}

int main(int argc, char **argv) {
	const char* output = "output.svg";
	auto [func, ground_truth] = function2d(argc,argv);
	std::size_t bins = 9;
	int plot_resolution = 100;
    unsigned long spp = 12;
    unsigned long spp_cv = 3;
	std::size_t seed = std::random_device()();
    double error_size_weight = 1.e-3;
    

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; } 
		else if (std::string(argv[i])=="-bins") { bins = atoi(argv[++i]); } 		
		else if (std::string(argv[i])=="-plot-resolution") { plot_resolution = atoi(argv[++i]); }
		else if (std::string(argv[i])=="-spp") { spp = atol(argv[++i]); } 		
		else if (std::string(argv[i])=="-spp-cv") { spp_cv = atol(argv[++i]); } 		
		else if (std::string(argv[i])=="-seed") { seed = atol(argv[++i]); }
        else if (std::string(argv[i])=="-error-size-weight") { error_size_weight=atof(argv[++i]); }
	}

	svg_cpp_plot::SVG svg;
	
	std::vector<double> gt(bins,0.0);
    double dx = 1.0/double(bins);
    for (std::size_t i=0; i<bins; ++i) gt[i]=double(bins)*ground_truth(dx*i,0,dx*(i+1),1);
	
    int hp=0;
    svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(hp,0)))
        .add(plot("Monte-Carlo",integrator_bins_stepper(stepper_bins_per_bin(stepper_monte_carlo_uniform(seed)),spp),func,range(0.0,0.0,1.0,1.0),bins,plot_resolution,gt));
//        .add(plot("Global Adaptive",stepper_bins_adaptive(nested(simpson,trapezoidal)),func,range(0.0,0.0,1.0,1.0),bins,max_samples,plot_resolution,gt));
//    svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(hp+=220,0)))
//        .add(plot("Global CV",stepper_bins_control_variate(control_variate_quadrature_adaptive(nested(simpson, trapezoidal),adaptive_iterations),seed),func,range(0.0,0.0,1.0,1.0),bins,max_samples,plot_resolution,gt));
//    svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(hp+=220,0)))
//        .add(plot("Region sampled CV",stepper_bins_adaptive_control_variates(nested(simpson, trapezoidal),adaptive_iterations,seed,seed+1),func,range(0.0,0.0,1.0,1.0),bins,max_samples,plot_resolution,gt));
	{
		auto& p = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(hp+=220,0)));
        p.add(plot("CV alpha=1",integrator_bins_stepper(stepper_bins_adaptive_stratified_control_variates(nested(simpson, trapezoidal),error_single_dimension_size(error_size_weight),(unsigned long)(bins*spp_cv/(3*2)),seed,seed+1), spp-spp_cv),func,range(0.0,0.0,1.0,1.0),bins,plot_resolution,gt));
		p.add(plot_adaptive(nested(simpson,trapezoidal),error_single_dimension_size(error_size_weight),func,range(0.0,0.0,1.0,1.0),(unsigned long)(bins*spp_cv/(3*2))));
	}
	
	{
		auto& p = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(hp+=220,0)));
        p.add(plot("CV optimized alpha",integrator_optimized_adaptive_stratified_control_variates(nested(simpson, trapezoidal),error_single_dimension_size(error_size_weight),(unsigned long)(bins*spp_cv/(3*2)),spp-spp_cv,seed),func,range(0.0,0.0,1.0,1.0),bins,plot_resolution,gt));
		p.add(plot_adaptive(nested(simpson,trapezoidal),error_single_dimension_size(error_size_weight),func,range(0.0,0.0,1.0,1.0),(unsigned long)(bins*spp_cv/(3*2))));
	}

	svg.viewBox(svg_cpp_plot::BoundingBox(-30,-240,hp+230,150));

	std::ofstream f(output);
	f << svg;
	
	return 0;
}



