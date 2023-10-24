#include "../quadrature/integrate.h"
#include "../functions/functions.h"
#include "../plot/integration2d.h"
#include <svg-cpp-plot/svg-cpp-plot.h>

#include <iostream>
#include <cmath>


void plot(const Function1D& func, int plot_resolution, const char* output) {
    auto [f,gt] = func;
    svg_cpp_plot::SVGPlot plt;
    plt.plot(svg_cpp_plot::arange(0.0,1.0,1.0/(plot_resolution)),f);
    plt.figsize({400,400});
    plt.savefig(output);
}

void plot(const Function2D& func, int plot_resolution, const char* output) {
    auto [f,gt] = func;
    float d = 1.0/float(plot_resolution);
    svg_cpp_plot::SVGPlot plt;
    plt.imshow(svg_cpp_plot::arange(0.5*d, 1.0, d), svg_cpp_plot::arange(0.5*d, 1.0, d), f).cmap("grayscale").interpolation("bilinear");
    plt.set_xticks({});
    plt.set_yticks({});
    plt.figsize({400,400});
    plt.savefig(output);
}

void plot(const Function3D& func, int plot_resolution, const char* output) {
    auto [f,gt] = func;
    float d = 1.0/float(plot_resolution);
    svg_cpp_plot::SVG svg;
    svg.viewBox({0,0,400,400});
    const unsigned int nzs = 5;
    float dz = 1.0f/float(nzs-1);
    float vmax = f(0,0,0); float vmin = f(0,0,0);
    for (float x = 0.5*d; x<1.0; x+=d)
        for (float y = 0.5*d; y<1.0; y+=d)
            for (float z = 0.5*d; z<1.0; z+=d) {
                float v = f(x,y,z);
                if (v>vmax) vmax = v;
                if (v<vmin) vmin = v;
            }
        
    for (unsigned int iz = 0; iz < nzs; ++iz) {
        float z = iz*dz;
        svg_cpp_plot::SVGPlot plt;
        plt.imshow(svg_cpp_plot::arange(0.5*d, 1.0, d), svg_cpp_plot::arange(0.5*d, 1.0, d), [f,z] (float x, float y) { return f(x,y,z); }).cmap("grayscale").interpolation("bilinear").vmax(vmax).vmin(vmin);
        plt.set_xticks({});
        plt.set_yticks({});
        plt.figsize({200,200});
        svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate({z*(400-200),z*(400-200)}))).add(plt.graph());
    }
    std::ofstream of(output); of<<svg;
}

void plot(const Function4D& func, int plot_resolution, const char* output) {
    auto [f,gt] = func;
    float d = 1.0/float(plot_resolution);
    svg_cpp_plot::SVG svg;
    svg.viewBox({0,0,400,400});
    const unsigned int nzs = 3;
    float dz = 1.0f/float(nzs-1);
    float dw = 1.0f/float(nzs-1);

    std::cout << "Loading vmax vmin" << std::endl;

    float vmax = f(0,0,0,0); float vmin = f(0,0,0,0);
    for (float x = 0.5*d; x<1.0; x+=d)
        for (float y = 0.5*d; y<1.0; y+=d)
            for (float z = 0.5*d; z<1.0; z+=d) {
                //for (float w = 0.5*d; w<1.0; w+=d) {
                float v = f(x,y,z,0);
                if (v>vmax) vmax = v;
                if (v<vmin) vmin = v;
            }

    for (unsigned int iz = 0; iz < nzs; ++iz)
    for (unsigned int iw = 0; iw < nzs; ++iw) {
        float z = iz*dz;
        float w = iw*dw;
        std::cout << "Ey " << z << "," << w << std::endl;
        svg_cpp_plot::SVGPlot plt;
        plt.imshow(svg_cpp_plot::arange(0.5*d, 1.0, d), svg_cpp_plot::arange(0.5*d, 1.0, d), [f,z,w] (float x, float y) { return f(x,y,z,w); }).cmap("grayscale").interpolation("bilinear").vmax(vmax).vmin(vmin);
        plt.set_xticks({});
        plt.set_yticks({});
        plt.figsize({133,133});
        svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate({iz*(133),iw*(133)}))).add(plt.graph());
    }
    std::ofstream of(output); of<<svg;
}

int main(int argc, char **argv) {
	const char* output = "output.svg";
	int iterations = 0;
	int plot_resolution = 100;

    auto func = function_from_commandline(argc,argv);

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; } 
		else if (std::string(argv[i])=="-iterations") { iterations = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-plot-resolution") { plot_resolution = atoi(argv[++i]); } 
	}
    
    std::visit([&] (auto&& f) { plot(f,plot_resolution,output); }, func);

	return 0;
}



