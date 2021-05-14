#include <memory>
#define MOTION_BLUR    //Usamos tres dimensiones (sólo podemos interpolar en la tercera y proyectar en las dos primeras)
#include "../../Method/division.h"
#include "../../Method/tupleMath.h"
#include "../../Method/multiArray.h"
#include "../../Method/errorScheduler.h"
#include <svg-cpp-plot/svg-cpp-plot.h>

#include <iostream>
#include <cmath>

typedef std::function<TupleMath(std::array<double, C_DIM_N>&, int, bool)> CALLBACK;

template<typename F,typename P>
class FunctionWrapper {
	F f;
	P& points;
public:
	FunctionWrapper(const F& f, P& points) : f(f), points(points) {}
	FunctionWrapper(F&& f, P& points) : f(std::forward<F>(f)), points(points) {}
	
	TupleMath operator()(const std::array<double, C_DIM_N>& pos, int n, bool b) const {
		double value = f(std::get<0>(pos),std::get<2>(pos)); //We ignore the 2nd image space dimension
		points.add_point(std::get<0>(pos),std::get<2>(pos));
		return TupleMath(value);
	}
};

template<typename F>
void plot(svg_cpp_plot::_2d::Group& group, const F& f, int maxIterations = 2, int initialImage = 1, int initialIllum = 1, int resolution = 200) {
	auto& graph_original = group.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(-250,0))).add(svg_cpp_plot::Graph2D({200,200},svg_cpp_plot::BoundingBox(0,0,1,1)));
	graph_original.area().add(svg_cpp_plot::_2d::function_2d(f,{0,0},{1,1},{resolution,resolution}));
	graph_original.border().stroke_width(1).stroke(svg_cpp_plot::black);

	auto& approximated_xt = group.add(svg_cpp_plot::_2d::group()).add(svg_cpp_plot::Graph2D({200,200},svg_cpp_plot::BoundingBox(0,0,1,1)));
	auto& graph_xt = approximated_xt.area().add(svg_cpp_plot::_2d::group());
	auto& regions_xt = approximated_xt.area().add(svg_cpp_plot::_2d::group());
	auto& approximated_xy = group.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,250))).add(svg_cpp_plot::Graph2D({200,200},svg_cpp_plot::BoundingBox(0,0,1,1)));
	auto& graph_xy = approximated_xy.area().add(svg_cpp_plot::_2d::group());
	auto& regions_xy = approximated_xy.area().add(svg_cpp_plot::_2d::group());


	approximated_xt.border().stroke_width(1).stroke(svg_cpp_plot::black);
	approximated_xy.border().stroke_width(1).stroke(svg_cpp_plot::black);
	
	auto& graph_error = group.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(250,0))).add(svg_cpp_plot::Graph2D({200,200},svg_cpp_plot::BoundingBox(0,0,1,1)));
	graph_error.border().stroke_width(1).stroke(svg_cpp_plot::black);

	auto points = svg_cpp_plot::_2d::points();
	FunctionWrapper<F,decltype(points)> wrapper(f,points);
	CALLBACK m_callback = wrapper;
	
	ErrorScheduler* scheduler = new ErrorScheduler("cacheDivisions/cacheW", 100, 1e-1, 25, 1.5, nullptr, C_DIM_N, C_DIM_L);

	/* Ajustar los parámetros de la división inicial */
	std::array<double, C_DIM_N> coords;
	std::array<double, C_DIM_N> gaps;
	coords.fill(0);
	gaps.fill(0.5);
	{
	    Division seedDivision(C_DIM_N, coords, gaps, nullptr, m_callback, false);
	    std::vector<std::shared_ptr<Division>> subDivisions = seedDivision.subdivideSeeds(initialImage, initialIllum, &(scheduler->pool));
	    for(auto ptr : subDivisions) scheduler->addDivisonError(ptr);
	}

    int iterations = 0;
    while(!scheduler->empty() && (maxIterations == -1 || iterations < maxIterations)) {
    	std::shared_ptr<Division> div = scheduler->getMaxError();
    	std::vector<std::shared_ptr<Division>> subDivisions = div->subdivide(2, true, &(scheduler->pool));
    	/* Por cada subdivision hijo */
    	for(std::shared_ptr<Division>& d : subDivisions) scheduler->addDivisonError(d);	
    	iterations++;
    }

    /* Si quedan divisiones por procesar */
    if (!scheduler->empty()) {
		double totalIntegral = 0;
    	for (auto & buff : scheduler->buffers) {
    		for(auto & d : buff) { 
    			// Calcula los coeficientes para la interpolación
    			auto coefs = d->interpolateCoefs(true);
    			// Calcula puntos 2d  (ignoramos una dimension, la 1)
			graph_xt.add(svg_cpp_plot::_2d::function_2d([&coefs,&d] (float x, float y) {
				return d->interpolateNDimensional(std::array<double,C_DIM_N>{0.5*(x - std::get<0>(d->coords))/std::get<0>(d->gaps), 0.0, 0.5*(y - std::get<2>(d->coords))/std::get<2>(d->gaps)})[0];
			}, {std::get<0>(d->coords),std::get<2>(d->coords)},
			   {std::get<0>(d->coords)+2.0*std::get<0>(d->gaps),std::get<2>(d->coords)+2.0*std::get<2>(d->gaps)},
			   {int(2*resolution*std::get<0>(d->gaps)),int(2*resolution*std::get<2>(d->gaps))}));
			
			graph_xy.add(svg_cpp_plot::_2d::function_2d([&coefs,&d] (float x, float y) {
				return d->interpolateNDimensional(std::array<double,C_DIM_N>{0.5*(x - std::get<0>(d->coords))/std::get<0>(d->gaps), 0.0, 0.5*(y - std::get<1>(d->coords))/std::get<1>(d->gaps)})[0];
			}, {std::get<0>(d->coords),std::get<1>(d->coords)},
			   {std::get<0>(d->coords)+2.0*std::get<0>(d->gaps),std::get<1>(d->coords)+2.0*std::get<1>(d->gaps)},
			   {int(2*resolution*std::get<0>(d->gaps)),int(2*resolution*std::get<1>(d->gaps))}));

			graph_error.area().add(svg_cpp_plot::_2d::function_2d([&f,&coefs,&d] (float x, float y) {
				return -1.0*(f(x,y) - d->interpolateNDimensional(std::array<double,C_DIM_N>{0.5*(x - std::get<0>(d->coords))/std::get<0>(d->gaps), 0.0, 0.5*(y - std::get<2>(d->coords))/std::get<2>(d->gaps)})[0]);
			}, svg_cpp_plot::_2d::color_map_red_blue(-1,1),
			   {std::get<0>(d->coords),std::get<2>(d->coords)},
			   {std::get<0>(d->coords)+2.0*std::get<0>(d->gaps),std::get<2>(d->coords)+2.0*std::get<2>(d->gaps)},
			   {int(2*resolution*std::get<0>(d->gaps)),int(2*resolution*std::get<2>(d->gaps))}));


			regions_xt.add(svg_cpp_plot::_2d::rect({std::get<0>(d->coords),std::get<2>(d->coords)},
			   {std::get<0>(d->coords)+2.0*std::get<0>(d->gaps),std::get<2>(d->coords)+2.0*std::get<2>(d->gaps)})).stroke(svg_cpp_plot::green).stroke_width(1);
			regions_xy.add(svg_cpp_plot::_2d::rect({std::get<0>(d->coords),std::get<1>(d->coords)},
			   {std::get<0>(d->coords)+2.0*std::get<0>(d->gaps),std::get<1>(d->coords)+2.0*std::get<1>(d->gaps)})).stroke(svg_cpp_plot::green).stroke_width(1);
			totalIntegral+=d->integral()[0];
    		}
    	}
		std::cout<<"TOTAL INTEGRAL = "<<totalIntegral<<std::endl;
    }


    approximated_xt.area().add(points).stroke_width(2).stroke(svg_cpp_plot::red);

    delete scheduler;
}



int main(int argc, char **argv) {
	QuadratureUtils::initStatic();	
	const char* output = "output.svg";
	int initial_subdivisions = 1;
	int iterations = 0;
	std::function<double(double,double)> func = [] (double x,double y) { return 0.5+0.4*(1-x)*std::sin(9*M_PI*x)*std::sin(3*M_PI*y); };
	int plot_resolution = 100;

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; } 
		else if (std::string(argv[i])=="-initial-subdivisions") { initial_subdivisions = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-iterations") { iterations = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-plot-resolution") { plot_resolution = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-function") { ++i;
			if (std::string(argv[i])=="step") {
				double atx = atof(argv[++i]);
				double aty = atof(argv[++i]);
				func = [atx,aty] (double x, double y) { return (((x/atx)+(y/aty))<1.0)?0.1:0.9; };
			}
			else if (std::string(argv[i])=="sin") {
				double freqx = atof(argv[++i]);
				double freqy = atof(argv[++i]);
				func = [freqx,freqy] (double x, double y) { return 0.5+0.4*std::sin(2*M_PI*x*freqx)*std::sin(2*M_PI*y*freqy); };
			}
			else if (std::string(argv[i])=="box") {
				double at1 = atof(argv[++i]); double at2 = atof(argv[++i]);
				func = [at1,at2] (double x, double y) { return (x<at1)?0.2:((y>at2)?0.7:0.8); };
			}
		}
	}

	svg_cpp_plot::SVG svg;
    	svg.viewBox(svg_cpp_plot::BoundingBox(-275,-25,525,525));
	plot(svg.add(svg_cpp_plot::_2d::group()),func,iterations, initial_subdivisions, initial_subdivisions, plot_resolution);
	std::ofstream f(output);
	f << svg;
	
	return 0;
}



