#include "../viltrum.h"
#include "../plot/integration2d.h"
#include <svg-cpp-plot/svg-cpp-plot.h>

#include <iostream>
#include <cmath>

using namespace viltrum;

int main(int argc, char **argv) {
	const char* output = "output.svg";
	int iterations = 0;
	std::function<double(double,double)> func = std::get<0>(function2d(argc,argv));
    auto f = [func] (const std::array<double,2>& x) { return func(x[0],x[1]); };
	int plot_resolution = 100;
    int width = 400;
    int spacing = 40;

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; } 
		else if (std::string(argv[i])=="-iterations") { iterations = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-plot-resolution") { plot_resolution = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-width") { width = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-spacing") { spacing = atoi(argv[++i]); } 
	}

    auto color_graph = svg_cpp_plot::rgb(0,0,0);
    auto color_map_function = svg_cpp_plot::_2d::color_map_grayscale();
//    auto color_map_function = svg_cpp_plot::_2d::color_map_heat(-1.0,1.0);
    auto color_map_error = svg_cpp_plot::_2d::color_map_red_blue(-1.0f,1.0f);
    auto color_cv  = svg_cpp_plot::rgb(0,0.6,0);
    float graph_width = 2;

    auto hp = 0;

	svg_cpp_plot::SVG svg;

    auto& graph_function = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,0))).add(svg_cpp_plot::Graph2D({width,width},svg_cpp_plot::BoundingBox(0,0,1,1)));
    graph_function.area().add(svg_cpp_plot::_2d::function_2d(func,color_map_function,{0.0f,0.0f},{1.0f,1.0f},{plot_resolution,plot_resolution}));
    graph_function.border().stroke_width(graph_width).stroke(color_graph);

	auto& graph_approx = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(hp+=(width+spacing),0))).add(svg_cpp_plot::Graph2D({width,width},svg_cpp_plot::BoundingBox(0,0,1,1)));
    graph_approx.area().add(plot_adaptive_approximation_2d(nested(simpson,trapezoidal),error_single_dimension_size(1.e-5f),f,range(0.0,0.0,1.0,1.0),iterations,color_map_function,plot_resolution));
    graph_approx.area().add(plot_adaptive_boundaries_2d(nested(simpson,trapezoidal),error_single_dimension_size(1.e-5f),f,range(0.0,0.0,1.0,1.0),iterations)).stroke_width(0.5*graph_width).stroke(color_cv);
    graph_approx.border().stroke_width(graph_width).stroke(color_graph);

	auto& graph_residual = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(hp+=(width+spacing),0))).add(svg_cpp_plot::Graph2D({width,width},svg_cpp_plot::BoundingBox(0,0,1,1)));
    graph_residual.area().add(plot_adaptive_error_2d(nested(simpson,trapezoidal),error_single_dimension_size(1.e-5f),f,range(0.0,0.0,1.0,1.0),iterations,color_map_error,plot_resolution));
    graph_residual.border().stroke_width(graph_width).stroke(color_graph);

    svg.viewBox(svg_cpp_plot::BoundingBox(-spacing,-spacing, hp+(width+spacing), width+spacing));
	std::ofstream file(output);
	file << svg;
	
	return 0;
}



