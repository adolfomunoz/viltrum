#include <memory>
//#define C_DIM_N 2 //Dos parece ser el mínimo o explota la interpolación polinomial. De todas formas ignoramos una de las dos dimensiones.
#include "../../Method/division.h"
#include "../../Method/tupleMath.h"
#include "../../Method/multiArray.h"
#include "../../Method/errorScheduler.h"
#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../quadrature/integrate.h"
#include "../functions/functions1d.h"

#include <iostream>
#include <cmath>

typedef std::function<TupleMath(std::array<double, C_DIM_N>&, int, bool)> CALLBACK;

template<typename F,typename P>
class FunctionWrapper {
	F f;
	P& points;
	mutable std::shared_ptr<unsigned long> nsamples;
public:
	FunctionWrapper(const F& f, P& points) : 
		f(f), points(points), nsamples(std::make_shared<unsigned long>(0)) {}
	FunctionWrapper(F&& f, P& points) : 
		f(std::forward<F>(f)), points(points), nsamples(std::make_shared<unsigned long>(0)) {}
	
	TupleMath operator()(const std::array<double, C_DIM_N>& pos, int n, bool b) const {
		double value = f(std::get<0>(pos));
		points.add_point(std::get<0>(pos),value);
		++(*nsamples);
		return TupleMath(value);
	}
	
	double operator()(const std::array<double,1>& pos) const {
		double value = f(std::get<0>(pos));
		points.add_point(std::get<0>(pos),value);
		++(*nsamples);
		return value;
	}
	
	unsigned long samples() const { return *nsamples; }
};

template<typename F>
void plot(svg_cpp_plot::Graph2D& graph, const F& f, int maxIterations = 2, int initialImage = 1, int initialIllum = 1) {
	graph.area().add(svg_cpp_plot::_2d::function(f,0,1,200)).fill(svg_cpp_plot::none).stroke_width(4).stroke(svg_cpp_plot::green);
	auto& points = graph.area().add(svg_cpp_plot::_2d::points()).stroke_width(10).stroke(svg_cpp_plot::black);
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
    			// Calcula puntos en el polinomio (ignoramos una dimension)
			graph.area().add(svg_cpp_plot::_2d::function([&coefs,&d] (float x) {
				return d->interpolatePoint(coefs, 0.5*(x - std::get<0>(d->coords))/std::get<0>(d->gaps), 0.0)[0];
			}, std::get<0>(d->coords),std::get<0>(d->coords)+2.0*std::get<0>(d->gaps),50)).stroke(svg_cpp_plot::red);

			graph.area().add(svg_cpp_plot::_2d::line({std::get<0>(d->coords),0},{std::get<0>(d->coords),1})).stroke_width(2).stroke(svg_cpp_plot::black).stroke_dasharray({2,2});
			graph.area().add(svg_cpp_plot::_2d::line({std::get<0>(d->coords)+2.0*std::get<0>(d->gaps),0},{std::get<0>(d->coords)+2.0*std::get<0>(d->gaps),1})).stroke_width(2).stroke(svg_cpp_plot::black).stroke_dasharray({2,2});
			totalIntegral+=d->integral()[0];
    		}
    	}
		
		std::stringstream label;
		label<<"Previous - "<<std::fixed<<std::setprecision(3)<<totalIntegral<<" ("<<wrapper.samples()<<" samples)";
		graph.add(svg_cpp_plot::_2d::text({400,-12},label.str())).font_size(38).text_anchor(svg_cpp_plot::text_anchor_middle);
    }

    delete scheduler;
}

template<typename F, typename N>
void plot_next(const char* name, svg_cpp_plot::Graph2D& graph, const F& f, const N& nested, int iterations = 2) {
	graph.area().add(svg_cpp_plot::_2d::function(f,0,1,200)).fill(svg_cpp_plot::none).stroke_width(4).stroke(svg_cpp_plot::green);
	auto& points = graph.area().add(svg_cpp_plot::_2d::points()).stroke_width(10).stroke(svg_cpp_plot::black);
	FunctionWrapper<F,decltype(points)> wrapper(f,points);
	
	auto regions = regions_adaptive_single_iterations(wrapper,nested,std::array{0.0},std::array{1.0},iterations);
	
	double totalIntegral(0);

	for (auto r : regions) {
		graph.area().add(svg_cpp_plot::_2d::function(
			[&r] (double x) { return r.approximation_at(std::array{x}); },
			std::get<0>(r.range())[0],std::get<1>(r.range())[0],50)).stroke(svg_cpp_plot::red);
		graph.area().add(svg_cpp_plot::_2d::line(
				{std::get<0>(r.range())[0],0},{std::get<0>(r.range())[0],1})).stroke_width(2).stroke(svg_cpp_plot::black).stroke_dasharray({2,2});
		graph.area().add(svg_cpp_plot::_2d::line(
				{std::get<1>(r.range())[0],0},{std::get<1>(r.range())[0],1})).stroke_width(2).stroke(svg_cpp_plot::black).stroke_dasharray({2,2});
		totalIntegral += r.integral();
	}
	std::stringstream label;
	label<<name<<" - "<<std::fixed<<std::setprecision(3)<<totalIntegral<<" ("<<wrapper.samples()<<" samples)";
	graph.add(svg_cpp_plot::_2d::text({400,-12},label.str())).font_size(38).text_anchor(svg_cpp_plot::text_anchor_middle);
}


int main(int argc, char **argv) {	
	const char* output = "output.svg";
	int initial_subdivisions = 1;
	int iterations = 0;
	auto func = std::get<0>(function1d(argc,argv));

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; } 
		else if (std::string(argv[i])=="-initial-subdivisions") { initial_subdivisions = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-iterations") { iterations = atoi(argv[++i]); } 
	}

	svg_cpp_plot::SVG svg;
    svg.viewBox(svg_cpp_plot::BoundingBox(-40,-40,840,1340));
	auto& graph_previous = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,0))).add(svg_cpp_plot::Graph2D({800,400},svg_cpp_plot::BoundingBox(0,0,1,1))); 
	auto& graph_boole_simpson = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,450))).add(svg_cpp_plot::Graph2D({800,400},svg_cpp_plot::BoundingBox(0,0,1,1)));	
	auto& graph_simpson_trapezoidal = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,900))).add(svg_cpp_plot::Graph2D({800,400},svg_cpp_plot::BoundingBox(0,0,1,1))); 

	plot(graph_previous,func, iterations, initial_subdivisions, initial_subdivisions);
	plot_next("Boole-Simpson",graph_boole_simpson,func,nested(boole,simpson), iterations);
	plot_next("Simpson-Trap.",graph_simpson_trapezoidal,func,nested(simpson, trapezoidal), iterations);
	std::ofstream f(output);
	f << svg;
	
	return 0;
}



