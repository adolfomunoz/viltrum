#include <memory>
#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../viltrum.h"
#include "../utils/cimg-wrapper.h"
#include <iostream>
#include <cmath>

int main(int argc, char **argv) {	
	const char* output = "motivation.svg";
//    std::size_t seed = std::random_device()();
    std::size_t seed = 7;
    float width = 654;
    float height = 480;
    long unsigned int spp = 3;
    int pixel_graphs = 4;
    int bins = 30;

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; }
		else if (std::string(argv[i])=="-seed") { seed = atol(argv[++i]); } 
		else if (std::string(argv[i])=="-width") { width = atof(argv[++i]); }
		else if (std::string(argv[i])=="-height") { height = atof(argv[++i]); }
        else if (std::string(argv[i])=="-spp") { spp = atoi(argv[++i]); }
        else if (std::string(argv[i])=="-pixel-graphs") { pixel_graphs = atoi(argv[++i]); }
        else if (std::string(argv[i])=="-bins") { bins = atol(argv[++i]); }
	}
    
    auto pixel = [] (float z) {return 1.0f - z*z; };
    auto f = [pixel] (float x, float y, float z) { return pixel(z); };

    auto color_integrand = svg_cpp_plot::rgb(1,0,0);
    auto color_samples = svg_cpp_plot::rgb(0,0,0);
    auto color_graph = svg_cpp_plot::rgb(0,0,0);
    float graph_width = height/100.0;
    float labelfontsize = height/10.0;
    

    {   
        auto pixel_integrator = viltrum::integrator_monte_carlo_uniform(spp,seed);
        auto image_integrator = viltrum::integrator_bins_per_bin(pixel_integrator);
        std::vector<std::vector<float>> data(bins,std::vector<float>(bins,0));
        viltrum::integrate_bins(image_integrator,data,viltrum::FunctionWrapper(f),viltrum::range(0.0f,0.0f,0.0f,1.0f,1.0f,1.0f));
        svg_cpp_plot::SVGPlot plt;
        plt.figsize({width,width}).xticks({}).yticks({});
        plt.set_xlabel("x").fontsize(labelfontsize);
        plt.set_ylabel("y").fontsize(labelfontsize);
        plt.imshow(data).cmap("grayscale").vmin(0).vmax(1);
        plt.savefig(std::string("montecarlo_")+output);
        for (int i = 0; i < pixel_graphs; ++i) {
            viltrum::FunctionWrapperProfile f(pixel);
            viltrum::integrate(pixel_integrator,f,viltrum::range(0.0f,1.0f));
            svg_cpp_plot::SVGPlot plt;
            plt.figsize({width,height}).xticks({}).yticks({}).axis({0,1,0,1});
            plt.plot(svg_cpp_plot::linspace(0.0f,1.0f,100),pixel).color(color_integrand).linewidth(graph_width);
            plt.scatter(f.params(0),f.values()).s(1.5*graph_width).c(color_samples);
            plt.set_xlabel("z").fontsize(labelfontsize);
            plt.savefig(std::string("montecarlo_")+std::to_string(i)+output);
        }
    }
    
    {   
        auto pixel_integrator = viltrum::integrator_quadrature(viltrum::simpson);
        auto image_integrator = viltrum::integrator_bins_per_bin(pixel_integrator);
        std::vector<std::vector<float>> data(bins,std::vector<float>(bins,0));
        viltrum::integrate_bins(image_integrator,data,viltrum::FunctionWrapper(f),viltrum::range(0.0f,0.0f,0.0f,1.0f,1.0f,1.0f));
        svg_cpp_plot::SVGPlot plt;
        plt.figsize({width,width}).xticks({}).yticks({});
        plt.set_xlabel("x").fontsize(labelfontsize);
        plt.set_ylabel("y").fontsize(labelfontsize);
        plt.imshow(data).cmap("grayscale").vmin(0).vmax(1);
        plt.savefig(std::string("trapezoidal_")+output);
        {
            viltrum::FunctionWrapperProfile f(pixel);
            viltrum::integrate(pixel_integrator,f,viltrum::range(0.0f,1.0f));
            svg_cpp_plot::SVGPlot plt;
            plt.figsize({width,height}).xticks({}).yticks({}).axis({0,1,0,1});
            plt.plot(svg_cpp_plot::linspace(0.0f,1.0f,100),pixel).color(color_integrand).linewidth(graph_width);
            plt.scatter(f.params(0),f.values()).s(1.5*graph_width).c(color_samples);
            plt.set_xlabel("z").fontsize(labelfontsize);
            plt.savefig(std::string("trapezoidal_0")+output);
            plt.plot(svg_cpp_plot::linspace(0.0f,1.0f,100),pixel).color("green").format("--").linewidth(1.1*graph_width);
            plt.savefig(std::string("trapezoidal_1")+output);
            plt.scatter(f.params(0),f.values()).s(1.5*graph_width).c(color_samples);
        }
    }

	return 0;
}



