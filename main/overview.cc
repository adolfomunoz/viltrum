#include <memory>
//#define C_DIM_N 2 //Dos parece ser el mínimo o explota la interpolación polinomial. De todas formas ignoramos una de las dos dimensiones.
#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../quadrature/integrate.h"
#include "../functions/functions1d.h"

#include <iostream>
#include <cmath>

template<typename R>
svg_cpp_plot::_2d::Group plot_regions(const std::vector<R>& regions) {
    auto sol = svg_cpp_plot::_2d::group();
    for (const auto& r : regions) 
        sol.add(svg_cpp_plot::_2d::function([r] (float x) { return r.approximation_at(std::array<float,1>{x}); }, r.range().min(0),r.range().max(0)));
    return sol;
}

template<typename F, typename R>
svg_cpp_plot::_2d::Group plot_residual(const F& f, const std::vector<R>& regions) {
    auto sol = svg_cpp_plot::_2d::group();
    for (const auto& r : regions) 
        sol.add(svg_cpp_plot::_2d::function([f,r] (float x) { return f(x) - r.approximation_at(std::array<float,1>{x}); }, r.range().min(0),r.range().max(0)));
    return sol;
}


template<typename R>
svg_cpp_plot::_2d::Group plot_boundaries(const std::vector<R>& regions) {
    auto sol = svg_cpp_plot::_2d::group();
    sol.add(svg_cpp_plot::_2d::line({0,0},{0,1}));
    for (const auto& r : regions) {
        sol.add(svg_cpp_plot::_2d::line({r.range().max(0),0},{r.range().max(0),1}));
    }
    return sol;
}


int main(int argc, char **argv) {	
	const char* output = "output.svg";
    float frequency = 3;
    float offset = 0.5;
    float scale = 0.6;
    std::size_t seed = std::random_device()();
	int iterations = 1;
    int iterations_end = 5;
    int width = 400;
    int spacing = 40;

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; }
        else if (std::string(argv[i])=="-frequency") { frequency = atof(argv[++i]); } 
		else if (std::string(argv[i])=="-iterations") { iterations = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-iterations-end") { iterations_end = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-seed") { seed = atol(argv[++i]); } 
		else if (std::string(argv[i])=="-width") { width = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-spacing") { spacing = atoi(argv[++i]); } 
	}

    auto color_integrand = svg_cpp_plot::rgb(1,0,0);
    auto color_pdf = svg_cpp_plot::rgb(0,0.5,1);
    auto color_cv  = svg_cpp_plot::rgb(0,0.6,0);
    auto color_graph = svg_cpp_plot::rgb(0,0,0);
    float graph_width = 2;

    auto pdf = [] (float x) { return 0.5*M_PI*std::cos(M_PI*(x - 0.5)); };
    auto cdf = [] (float x) { return 0.5*std::sin(M_PI*(x - 0.5)) + 0.5; };
    auto integrand_primary = [&] (float x) { 
        return scale*(siv::PerlinNoise(std::uint32_t(seed)).noise1D(x*frequency) +
                  0.5*siv::PerlinNoise(std::uint32_t(seed+1)).noise1D(10.0f*x*frequency))  + offset; };
    auto integrand = [&] (float x) { return pdf(x)*integrand_primary(cdf(x)); };
    auto plottable_pdf = [&] (float x) { return pdf(x)/(0.7*M_PI); };
    auto adaptable_integrand = [&] (const std::array<float,1>& x) { return integrand_primary(x[0]); };

    auto hp = 0;

	svg_cpp_plot::SVG svg;

    auto& graph_previous = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,0))).add(svg_cpp_plot::Graph2D({width,width},svg_cpp_plot::BoundingBox(0,0,1,1)));
    graph_previous.area().add(svg_cpp_plot::_2d::function(plottable_pdf,0,1)).stroke_width(graph_width).stroke(color_pdf); 
    graph_previous.area().add(svg_cpp_plot::_2d::function(integrand,0,1)).stroke_width(graph_width).stroke(color_integrand); 
    graph_previous.area().add(svg_cpp_plot::_2d::line({0,0},{0,1})).stroke_width(graph_width).stroke(color_graph); 
    graph_previous.area().add(svg_cpp_plot::_2d::line({0,0},{1,0})).stroke_width(graph_width).stroke(color_graph);

	auto& graph_primary = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(hp+=(width+spacing),0))).add(svg_cpp_plot::Graph2D({width,width},svg_cpp_plot::BoundingBox(0,0,1,1)));
    graph_primary.area().add(svg_cpp_plot::_2d::function(integrand_primary,0,1)).stroke_width(graph_width).stroke(color_integrand); 
    graph_primary.area().add(svg_cpp_plot::_2d::line({0,0},{0,1})).stroke_width(graph_width).stroke(color_graph); 
    graph_primary.area().add(svg_cpp_plot::_2d::line({0,0},{1,0})).stroke_width(graph_width).stroke(color_graph);

    auto stepper = stepper_adaptive(nested(simpson,trapezoidal));
    auto regions = stepper.init(adaptable_integrand,range(0.0f,1.0f));
    for (int i = 0; i < iterations; ++i) stepper.step(adaptable_integrand, range(0.0f,1.0f), regions);

    for (int i = 0; i<2; ++i) {
        auto& graph_cv  = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(hp+=(width+spacing),0))).add(svg_cpp_plot::Graph2D({width,width},svg_cpp_plot::BoundingBox(0,0,1,1)));
        graph_cv.area().add(plot_regions(regions)).stroke_width(graph_width).stroke(color_cv);
        graph_cv.area().add(svg_cpp_plot::_2d::function(integrand_primary,0,1)).stroke_width(graph_width).stroke(color_integrand); 
        graph_cv.area().add(plot_boundaries(regions)).stroke_width(0.5*graph_width).stroke(color_cv);
        graph_cv.area().add(svg_cpp_plot::_2d::line({0,0},{0,1})).stroke_width(graph_width).stroke(color_graph); 
        graph_cv.area().add(svg_cpp_plot::_2d::line({0,0},{1,0})).stroke_width(graph_width).stroke(color_graph);
        stepper.step(adaptable_integrand, range(0.0f,1.0f), regions);
    }

    for (int i = 0; i < (iterations_end - iterations - 1); ++i) stepper.step(adaptable_integrand, range(0.0f,1.0f), regions);

    auto& graph_cv  = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(hp+=(width+spacing),0))).add(svg_cpp_plot::Graph2D({width,width},svg_cpp_plot::BoundingBox(0,0,1,1)));
    graph_cv.area().add(plot_regions(regions)).stroke_width(graph_width).stroke(color_cv);
    graph_cv.area().add(svg_cpp_plot::_2d::function(integrand_primary,0,1)).stroke_width(graph_width).stroke(color_integrand); 
    graph_cv.area().add(svg_cpp_plot::_2d::line({0,0},{0,1})).stroke_width(graph_width).stroke(color_graph); 
    graph_cv.area().add(svg_cpp_plot::_2d::line({0,0},{1,0})).stroke_width(graph_width).stroke(color_graph);


   auto& graph_residual  = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(hp+=(width+spacing),0))).add(svg_cpp_plot::Graph2D({width,width},svg_cpp_plot::BoundingBox(0,-0.5,1,0.5)));
    graph_residual.area().add(svg_cpp_plot::_2d::line({0,0},{1,0})).stroke_width(graph_width).stroke(color_cv);
    graph_residual.area().add(plot_residual(integrand_primary,regions)).stroke_width(graph_width).stroke(color_integrand); 
    graph_residual.area().add(svg_cpp_plot::_2d::line({0,-0.5},{0,0.5})).stroke_width(graph_width).stroke(color_graph);
    
    svg.viewBox(svg_cpp_plot::BoundingBox(-spacing,-spacing, hp+(width+spacing), width+spacing));
	std::ofstream f(output);
	f << svg;
	
	return 0;
}



