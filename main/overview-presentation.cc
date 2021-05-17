#include <memory>
#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../viltrum.h"
#include "../functions/functions1d.h"

#include <iostream>
#include <cmath>

int main(int argc, char **argv) {	
	const char* output = "overview.svg";
    float frequency = 1;
    float offset = 0.5;
    float scale = 0.6;
//    std::size_t seed = std::random_device()();
    std::size_t seed = 2;
	int iterations = 9;
    float width = 1280;
    float height = 720;
    int bins = 10;
    int plot_samples = 1000;

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; }
        else if (std::string(argv[i])=="-frequency") { frequency = atof(argv[++i]); } 
		else if (std::string(argv[i])=="-iterations") { iterations = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-seed") { seed = atol(argv[++i]); } 
		else if (std::string(argv[i])=="-width") { width = atof(argv[++i]); }
		else if (std::string(argv[i])=="-height") { height = atof(argv[++i]); }
        else if (std::string(argv[i])=="-bins") { bins = atoi(argv[++i]); }
        else if (std::string(argv[i])=="-plot-samples") { plot_samples = atoi(argv[++i]); }
	}

    auto color_integrand = svg_cpp_plot::rgb(1,0,0);
    auto color_pdf = svg_cpp_plot::rgb(0,0.5,1);
    auto color_cv  = svg_cpp_plot::rgb(0,0.6,0);
    auto color_graph = svg_cpp_plot::rgb(0,0,0);
    float graph_width = height/100.0;
    float tickfontsize = height/20.0;

    auto pdf = [] (float x) { return 0.5*M_PI*std::cos(M_PI*(x - 0.5)); };
    auto cdf = [] (float x) { return 0.5*std::sin(M_PI*(x - 0.5)) + 0.5; };
    auto integrand_primary = [&] (float x) { 
        return scale*(siv::PerlinNoise(std::uint32_t(seed)).noise1D(x*frequency)  + offset; };
    auto integrand = [&] (float x) { return pdf(x)*integrand_primary(cdf(x)); };
    auto plottable_pdf = [&] (float x) { return pdf(x)/(0.7*M_PI); };
    auto adaptable_integrand = [&] (const std::array<float,1>& x) { return integrand_primary(x[0]); };

    {   //Integrand
        svg_cpp_plot::SVGPlot plt;
        plt.figsize({width,height});
        plt.xticks({0,1}).set_xticklabels({"a","b"}).fontsize(tickfontsize);
        plt.set_yticks({0}).fontsize(tickfontsize);
        plt.plot(svg_cpp_plot::linspace(0,1,plot_samples),integrand).linewidth(graph_width).color(color_integrand);
        plt.savefig(std::string("integrand_")+output);
        //Integrand+pdf
        plt.plot(svg_cpp_plot::linspace(0,1,plot_samples),plottable_pdf).linewidth(graph_width).color(color_pdf);
        plt.savefig(std::string("integrand_pdf_")+output);
    }
    
    {   //Primary space
        svg_cpp_plot::SVGPlot plt;
        plt.figsize({width,height}).axis({0,1,0,1});
        plt.set_xticks({0,1}).fontsize(tickfontsize);
        plt.set_yticks({0}).fontsize(tickfontsize);
        plt.plot(svg_cpp_plot::linspace(0,1,plot_samples),integrand_primary).linewidth(graph_width).color(color_integrand);
        plt.savefig(std::string("primary_")+output);
    }
    auto stepper = viltrum::stepper_adaptive(viltrum::nested(viltrum::simpson,viltrum::trapezoidal));
    auto regions = stepper.init(adaptable_integrand,viltrum::range(0.0f,1.0f));
    for (int i = 0; i < iterations; ++i) { //Iterations
        svg_cpp_plot::SVGPlot plt;
        plt.figsize({width,height}).axis({0,1,0,1});
        plt.set_xticks({0,1}).fontsize(tickfontsize);
        plt.set_yticks({0}).fontsize(tickfontsize);
        plt.plot(svg_cpp_plot::linspace(0,1,plot_samples),integrand_primary).linewidth(graph_width).color(color_integrand);
        for (const auto& r : regions) {
            plt.plot(svg_cpp_plot::linspace(r.range().min(0),r.range().max(0),float(plot_samples)/(r.range().max(0)-r.range().min(0))),
                [r] (float x) { return r.approximation_at(std::array<float,1>{x}); })
                    .linewidth(graph_width).color(color_cv);
            plt.plot({r.range().min(0),r.range().min(0)},{0,1}).linewidth(0.5*graph_width).color(color_cv);
            plt.plot({r.range().max(0),r.range().max(0)},{0,1}).linewidth(0.5*graph_width).color(color_cv); 
        }
        plt.savefig(std::string("iteration")+std::to_string(i)+"_"+output);
        if (i < (iterations-1)) stepper.step(adaptable_integrand, viltrum::range(0.0f,1.0f), regions);
    }
        
    {
        svg_cpp_plot::SVGPlot plt;
        plt.figsize({width,height}).axis({0,1,-0.5,0.5});
        plt.set_xticks({0,1}).fontsize(tickfontsize);
        plt.set_yticks({0}).fontsize(tickfontsize);
        plt.plot({0,1},{0,0}).linewidth(0.5*graph_width).color(color_graph);
        for (const auto& r : regions) {
            plt.plot(
                svg_cpp_plot::linspace(r.range().min(0),r.range().max(0),float(plot_samples)/(r.range().max(0)-r.range().min(0))),
                [r,integrand_primary] (float x) { return integrand_primary(x) - r.approximation_at(std::array<float,1>{x}); })
                    .linewidth(graph_width).color(color_integrand);
        }
        
        plt.savefig(std::string("residual_")+output);
    }

	return 0;
}



