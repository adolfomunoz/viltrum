#include <iostream>
#include "../multiarray/fold.h"
#include "../multiarray/fill.h"
#include "../multiarray/multiarray.h"
#include "../quadrature/rules.h"
#include <cmath>
#include <random>
#include <svg-cpp-plot/svg-cpp-plot.h>

template<typename F,typename P>
class FunctionWrapper {
	F f;
	P& points;
public:
	FunctionWrapper(const F& f, P& points) : f(f), points(points) {}
	FunctionWrapper(F&& f, P& points) : f(std::forward<F>(f)), points(points) {}
	
	double operator()(const std::array<double, 1>& pos) const {
		double value = f(std::get<0>(pos));
		points.add_point(std::get<0>(pos),value);
		return value;
	}
};

template<typename F, typename Q>
void test(const char* name, svg_cpp_plot::Graph2D& graph, const F& f, const Q& quadrature) {
	graph.border().stroke_width(1).stroke(svg_cpp_plot::black);
	graph.xticks(2).stroke_width(1).stroke(svg_cpp_plot::black);
	graph.xlabels(2);


	graph.area().add(svg_cpp_plot::_2d::function(f,0,1,100)).fill(svg_cpp_plot::none).stroke_width(4).stroke(svg_cpp_plot::green);
	auto& points = graph.area().add(svg_cpp_plot::_2d::points()).stroke_width(10).stroke(svg_cpp_plot::black);
	FunctionWrapper<F,decltype(points)> wrapper(f,points);

	multiarray<double,Q::samples,1> ma;
	ma.fill(wrapper);
//	std::cout<<ma.fold(array_to_string(" ")).value()<<std::endl;
/*	double mc = 0;
	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with random device
    	std::uniform_real_distribution<double> dis(0.0, 1.0);
    	for (int n=0; n<10000; ++n) mc+=f(dis(gen));
	mc/=10000.0;*/
	
	std::stringstream label;
	label<<name<<" ("<<std::fixed<<std::setprecision(3)<<ma.fold_all(quadrature)<<")";
	graph.add(svg_cpp_plot::_2d::text({200,-12},label.str())).font_size(38).text_anchor(svg_cpp_plot::text_anchor_middle);


//	std::cout<<name<<"\t Integral = "<<fold_all(ma,quadrature)<<" - "<<mc<<std::endl;
	graph.area().add(svg_cpp_plot::_2d::function([ma,quadrature] (double t) -> double { return ma.fold([&] (const auto& v) { return quadrature.at(t,v); }); },0,1,100)).stroke_width(2).stroke(svg_cpp_plot::red);


}


int main(int argc, char **argv) {
	std::function<double(double)> f = [] (double x) { return 0.5*std::cos(M_PI*x/2.0) + 0.5; };
//	auto f = [] (double x) { return x; };

	const char* output = "output.svg";
	for (int i = 0; i < argc; ++i) {
		if ( (std::string(argv[i]) == "-output") && (i<(argc-1)) )
			output = argv[++i];
		else if ( (std::string(argv[i]) == "-function") && (i<(argc-2)) ) {
			std::string function = argv[++i];
			double param = atof(argv[++i]);
			if (function == "cos") 
				f = [param] (double x) { return 0.5*std::cos(2*M_PI*x/param)+0.5; };
			else if (function == "step") {
				f = [param] (double x) { return (x<param)?0.2:0.8; };
			}
		}
	}

	svg_cpp_plot::SVG svg;
    	svg.viewBox(svg_cpp_plot::BoundingBox(-40,-40,1320,960));
	
	test("Trapezoidal  ",
		svg.add(svg_cpp_plot::_2d::group())
			.add(svg_cpp_plot::Graph2D({400,400},svg_cpp_plot::BoundingBox(0,0,1,1))),
			f,trapezoidal);
	test("Simpson      ",
		svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate({440,0})))
			.add(svg_cpp_plot::Graph2D({400,400},svg_cpp_plot::BoundingBox(0,0,1,1))),
			f,simpson);
	test("Boole        ",
		svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate({880,0})))
			.add(svg_cpp_plot::Graph2D({400,400},svg_cpp_plot::BoundingBox(0,0,1,1))),
			f,boole);
	test("Trapezoidal x2",
		svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate({0,480})))
			.add(svg_cpp_plot::Graph2D({400,400},svg_cpp_plot::BoundingBox(0,0,1,1))),
			f,steps<2>(trapezoidal));
	test("Trapezoidal x4",
		svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate({440,480})))
			.add(svg_cpp_plot::Graph2D({400,400},svg_cpp_plot::BoundingBox(0,0,1,1))),
			f,steps<4>(trapezoidal));
	test("Simpson x2    ",
		svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate({880,480})))
			.add(svg_cpp_plot::Graph2D({400,400},svg_cpp_plot::BoundingBox(0,0,1,1))),
			f,steps<2>(simpson));	
	std::ofstream of(output);
	of << svg;
	return 0;
}



