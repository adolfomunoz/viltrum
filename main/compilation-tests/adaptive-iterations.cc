#include "../../viltrum.h"
#include <iostream>
#include <iomanip>
#include <string>

using namespace viltrum;


template<typename Integrator, typename Function>
void test(const std::string& name, const Integrator& integrator, const Function& f, const std::size_t bins) {
    LoggerProgress logger(name);
    std::vector<float> sol(bins,0.0f);
    integrate(integrator,sol,f,range_primary<2>(),logger);
    for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
    std::cout<<std::endl;
}

int main(int argc, char** argv) {
    std::size_t bins = 24;
    unsigned long iterations = 100000;

    for (int i=1;i<argc;++i) {
        std::string arg(argv[i]);
        if (arg=="-bins") {
            if (i+1<argc) {
                bins = std::stoul(argv[++i]);
            } else {
                std::cerr<<"Missing argument for -bins"<<std::endl;
                return 1;
            }
        } else if (arg=="-iterations") {
            if (i+1<argc) {
                iterations = std::stoul(argv[++i]);
            } else {
                std::cerr<<"Missing argument for -iterations"<<std::endl;
                return 1;
            }
        } else {
            std::cerr<<"Unknown argument "<<arg<<std::endl;
            return 1;
        }
    }

    auto f =[] (const std::array<float,2>& x) -> double {
       if ((x[0]+x[1])<1) return 1.0;        
       else return 0.0;
    };

    auto error_heuristic = error_heuristic_default(error_metric_absolute());

    std::cout<<"Integration of a 2D function with adaptive iterations, bins="<<bins<<", iterations="<<iterations<<std::endl;

    test("Sequential                       ",
        integrator_adaptive_iterations(nested(simpson,trapezoidal),error_heuristic, iterations),
        f,bins);
    
    test("Parallel - multiqueue  | regions ",
        integrator_adaptive_iterations_parallel(nested(simpson,trapezoidal),error_heuristic, heap::multiqueue(), iterations),
        f,bins);

    test("Parallel - multiqueue  | bins    ",
        integrator_adaptive_iterations_parallel_bins(nested(simpson,trapezoidal),error_heuristic, heap::multiqueue(), iterations),
        f,bins);

    test("Parallel - lockedqueue | regions ",
        integrator_adaptive_iterations_parallel(nested(simpson,trapezoidal),error_heuristic, heap::lockedpriorityqueue(), iterations),
        f,bins);

    test("Parallel - lockedqueue | bins    ",
        integrator_adaptive_iterations_parallel_bins(nested(simpson,trapezoidal),error_heuristic, heap::lockedpriorityqueue(), iterations),
        f,bins);

} 