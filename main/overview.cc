#include <memory>
//#define C_DIM_N 2 //Dos parece ser el mínimo o explota la interpolación polinomial. De todas formas ignoramos una de las dos dimensiones.
#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../quadrature/integrate.h"
#include "../functions/functions1d.h"

#include <iostream>
#include <cmath>

/*
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
*/


int main(int argc, char **argv) {	
	const char* output = "output.svg";
    float frequency = 3;
    float offset = 0.5;
    float scale = 0.6;
    std::size_t seed = std::random_device()();
	int iterations = 1;
    int iterations_end = 5;
    float width = 200;
    int plot_samples = 100;
//    int spacing = 40;

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; }
        else if (std::string(argv[i])=="-frequency") { frequency = atof(argv[++i]); } 
		else if (std::string(argv[i])=="-iterations") { iterations = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-iterations-end") { iterations_end = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-seed") { seed = atol(argv[++i]); } 
		else if (std::string(argv[i])=="-width") { width = atof(argv[++i]); }
        else if (std::string(argv[i])=="-plot-samples") { plot_samples = atoi(argv[++i]); }
//		else if (std::string(argv[i])=="-spacing") { spacing = atoi(argv[++i]); } 
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
                  0.35*siv::PerlinNoise(std::uint32_t(seed+1)).noise1D(5.4f*x*frequency))  + offset; };
    auto integrand = [&] (float x) { return pdf(x)*integrand_primary(cdf(x)); };
    auto plottable_pdf = [&] (float x) { return pdf(x)/(0.7*M_PI); };
    auto adaptable_integrand = [&] (const std::array<float,1>& x) { return integrand_primary(x[0]); };

	svg_cpp_plot::SVG svg;

    svg_cpp_plot::SVGPlot plt;
    for (int i = 0; i<6; ++i) plt.subplot(1,6,i).figsize({width,width}).xticks({0,1}).yticks({0}).axis({0,1,0,1});
    for (int i = 1; i<5; ++i) plt.subplot(1,6,i).axis({0,1,0,1});
    plt.subplots_adjust().wspace(0.02);
    plt.subplot(1,6,0).plot(svg_cpp_plot::linspace(0,1,plot_samples),plottable_pdf).linewidth(graph_width).color(color_pdf);
    plt.subplot(1,6,0).plot(svg_cpp_plot::linspace(0,1,plot_samples),integrand).linewidth(graph_width).color(color_integrand);
    plt.subplot(1,6,0).xticklabels({"a","b"});
    plt.subplot(1,6,1).plot(svg_cpp_plot::linspace(0,1,plot_samples),integrand_primary).linewidth(graph_width).color(color_integrand);

    auto stepper = stepper_adaptive(nested(simpson,trapezoidal));
    auto regions = stepper.init(adaptable_integrand,range(0.0f,1.0f));
    for (int i = 0; i < iterations; ++i) stepper.step(adaptable_integrand, range(0.0f,1.0f), regions);

    for (int i = 0; i<=2; ++i) {
        for (const auto& r : regions) {
            plt.subplot(1,6,2+i).plot(
                svg_cpp_plot::linspace(r.range().min(0),r.range().max(0),float(plot_samples)/(r.range().max(0)-r.range().min(0))),
                [r] (float x) { return r.approximation_at(std::array<float,1>{x}); })
                    .linewidth(graph_width).color(color_cv);
            plt.subplot(1,6,2+i).plot({r.range().min(0),r.range().min(0)},{0,1})
                    .linewidth(0.5*graph_width).color(color_cv);
            plt.subplot(1,6,2+i).plot({r.range().max(0),r.range().max(0)},{0,1})
                    .linewidth(0.5*graph_width).color(color_cv);
        }

        plt.subplot(1,6,2+i).plot(svg_cpp_plot::linspace(0,1,100),integrand_primary).linewidth(graph_width).color(color_integrand);
        if (i<2) stepper.step(adaptable_integrand, range(0.0f,1.0f), regions);
        if (i==1) for(int j=i;j<(iterations_end-iterations)-1;++j) stepper.step(adaptable_integrand, range(0.0f,1.0f), regions);
    }
    
    for (const auto& r : regions) {
        plt.subplot(1,6,5).plot(
            svg_cpp_plot::linspace(r.range().min(0),r.range().max(0),float(plot_samples)/(r.range().max(0)-r.range().min(0))),
            [r,integrand_primary] (float x) { return integrand_primary(x) - r.approximation_at(std::array<float,1>{x}); })
                    .linewidth(graph_width).color(color_integrand);
    }
    plt.subplot(1,6,5).axis({0,1,-0.5,0.5}).plot({0,1},{0,0}).linewidth(0.5*graph_width).color(color_graph);

    plt.savefig(output);   
	
	return 0;
}



