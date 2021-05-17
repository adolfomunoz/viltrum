#include <memory>
#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../viltrum.h"
#include "../functions/functions1d.h"

#include <iostream>
#include <cmath>

int main(int argc, char **argv) {	
	const char* output = "overview.svg";
    float frequency = 3;
    float offset = 0.5;
    float scale = 0.6;
    std::size_t seed = std::random_device()();
	int iterations = 9;
    float width = 1280;
    float height = 720;
    int plot_samples = 100;
//    int spacing = 40;

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; }
        else if (std::string(argv[i])=="-frequency") { frequency = atof(argv[++i]); } 
		else if (std::string(argv[i])=="-iterations") { iterations = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-seed") { seed = atol(argv[++i]); } 
		else if (std::string(argv[i])=="-width") { width = atof(argv[++i]); }
		else if (std::string(argv[i])=="-height") { height = atof(argv[++i]); }
        else if (std::string(argv[i])=="-plot-samples") { plot_samples = atoi(argv[++i]); }
	}

    auto color_integrand = svg_cpp_plot::rgb(1,0,0);
    auto color_pdf = svg_cpp_plot::rgb(0,0.5,1);
    auto color_cv  = svg_cpp_plot::rgb(0,0.6,0);
    auto color_graph = svg_cpp_plot::rgb(0,0,0);
    float graph_width = height/100.0;
    float tickfontsize = height/20.0;

    auto integrand_primary = [&] (float x) { return 1.0f/(1.0f + 25.0f*x*x); };

    {   //Simpson
        svg_cpp_plot::SVGPlot plt;
        viltrum::FunctionWrapperProfile f(integrand_primary);
        auto stepper = viltrum::stepper_adaptive(viltrum::nested(viltrum::simpson,viltrum::trapezoidal));
        auto regions = stepper.init(f,viltrum::range(-1.0f,1.0f));
        plt.figsize({width,height}).axis({-1,1,-1,1});
        plt.xticks({});
//        plt.set_xticks({-1,1}).fontsize(tickfontsize);
        plt.set_yticks({0}).fontsize(tickfontsize);
        plt.plot(svg_cpp_plot::linspace(-1,1,plot_samples),integrand_primary).linewidth(graph_width).color(color_integrand);
        for (const auto& r : regions) {
            plt.plot(svg_cpp_plot::linspace(r.range().min(0),r.range().max(0),float(plot_samples)/(r.range().max(0)-r.range().min(0))),
                [r] (float x) { return r.approximation_at(std::array<float,1>{x}); })
                    .linewidth(graph_width).color(color_cv);
        }
        plt.scatter(f.params(0),f.values()).c(color_integrand).edgecolors(color_cv).s(2*graph_width).linewidths(0.5*graph_width);
        plt.savefig(std::string("simpson_")+output);
    }
    {   //Bool
        svg_cpp_plot::SVGPlot plt;
        viltrum::FunctionWrapperProfile f(integrand_primary);
        auto stepper = viltrum::stepper_adaptive(viltrum::nested(viltrum::boole,viltrum::simpson));
        auto regions = stepper.init(f,viltrum::range(-1.0f,1.0f));
        plt.figsize({width,height}).axis({-1,1,-1,1});
        plt.xticks({});
//        plt.set_xticks({-1,1}).fontsize(tickfontsize);
        plt.set_yticks({0}).fontsize(tickfontsize);
        plt.plot(svg_cpp_plot::linspace(-1,1,plot_samples),integrand_primary).linewidth(graph_width).color(color_integrand);
        for (const auto& r : regions) {
            plt.plot(svg_cpp_plot::linspace(r.range().min(0),r.range().max(0),float(plot_samples)/(r.range().max(0)-r.range().min(0))),
                [r] (float x) { return r.approximation_at(std::array<float,1>{x}); })
                    .linewidth(graph_width).color(color_cv);
        }
        plt.scatter(f.params(0),f.values()).c(color_integrand).edgecolors(color_cv).s(2*graph_width).linewidths(0.5*graph_width);
        plt.savefig(std::string("boole_")+output);
    }

	return 0;
}



