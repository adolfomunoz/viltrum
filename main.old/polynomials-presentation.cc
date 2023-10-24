#include <memory>
#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../viltrum.h"
#include "../functions/functions1d.h"

#include <iostream>
#include <cmath>

int main(int argc, char **argv) {	
	const char* output = "polynomials.svg";
	std::size_t bins = 9;
    float width = 1280;
    float height = 720;
    int plot_samples = 100;

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; }
		else if (std::string(argv[i])=="-bins") { bins = atol(argv[++i]); } 
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

    auto integrand_runge = [&] (float x) { return 1.0f/(1.0f + 25.0f*x*x); };
    auto integrand_smooth =[&] (float x) { return x/(1.0f+std::exp(x)); };


    {   //Bins
        svg_cpp_plot::SVGPlot plt;
        std::vector<float> solv(bins);
        auto sol = viltrum::adaptor(solv);
        viltrum::FunctionWrapperProfile f(integrand_smooth);
        auto stepper = viltrum::stepper_bins_adaptive(viltrum::nested(viltrum::simpson,viltrum::trapezoidal));
        auto regions = stepper.init(sol.resolution(),f,viltrum::range(0.0f,1.0f));
        plt.figsize({width,height});
        plt.xticks({});
//        plt.set_xticks({-1,1}).fontsize(tickfontsize);
        plt.set_yticks({0}).fontsize(tickfontsize);
        plt.plot(svg_cpp_plot::linspace(0,1,plot_samples),integrand_smooth).linewidth(graph_width).color(color_integrand);
        for (const auto& r : regions) {
            plt.plot(svg_cpp_plot::linspace(r.range().min(0),r.range().max(0),float(plot_samples)/(r.range().max(0)-r.range().min(0))),
                [r] (float x) { return r.approximation_at(std::array<float,1>{x}); })
                    .linewidth(graph_width).color(color_cv);
        }
        plt.scatter(f.params(0),f.values()).c(color_integrand).edgecolors(color_cv).s(2*graph_width).linewidths(0.5*graph_width);
        plt.savefig(std::string("approx_")+output);
        
        svg_cpp_plot::SVGPlot plti;
        plti.figsize({width,height}).xticks({}).set_yticks({0}).fontsize(tickfontsize);
        stepper.integral(sol, sol.resolution(), f, viltrum::range(0.0f,1.0f), regions); 
        plti.bar(svg_cpp_plot::arange(bins),solv).width({1}).color(svg_cpp_plot::rgb(0.75,0.75,0));
        plti.savefig(std::string("bins_")+output);
    }
    
    {   //Simpson
        svg_cpp_plot::SVGPlot plt;
        viltrum::FunctionWrapperProfile f(integrand_runge);
        auto stepper = viltrum::stepper_adaptive(viltrum::nested(viltrum::simpson,viltrum::trapezoidal));
        auto regions = stepper.init(f,viltrum::range(-1.0f,1.0f));
        stepper.step(f,viltrum::range(-1.0f,1.0f),regions);
        plt.figsize({width,height}).axis({-1,1,-1,1});
        plt.xticks({});
//        plt.set_xticks({-1,1}).fontsize(tickfontsize);
        plt.set_yticks({0}).fontsize(tickfontsize);
        plt.plot(svg_cpp_plot::linspace(-1,1,plot_samples),integrand_runge).linewidth(graph_width).color(color_integrand);
        plt.plot({-1,1},{0,0}).linewidth(0.5*graph_width).color(color_graph);
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
        viltrum::FunctionWrapperProfile f(integrand_runge);
        auto stepper = viltrum::stepper_adaptive(viltrum::nested(viltrum::boole,viltrum::simpson));
        auto regions = stepper.init(f,viltrum::range(-1.0f,1.0f));
        plt.figsize({width,height}).axis({-1,1,-1,1});
        plt.xticks({});
//        plt.set_xticks({-1,1}).fontsize(tickfontsize);
        plt.set_yticks({0}).fontsize(tickfontsize);
        plt.plot(svg_cpp_plot::linspace(-1,1,plot_samples),integrand_runge).linewidth(graph_width).color(color_integrand);
        plt.plot({-1,1},{0,0}).linewidth(0.5*graph_width).color(color_graph);
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



