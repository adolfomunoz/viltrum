#pragma once

#include "functions1d.h"
#include "functions2d.h"
#include "functions3d.h"
#include "functions4d.h"

#include <variant>

using Function = std::variant<Function1D, Function2D, Function3D, Function4D>;

std::size_t function_dimensions(const Function1D& f) {
    return 1;
}
std::size_t function_dimensions(const Function2D& f) {
    return 2;
}
std::size_t function_dimensions(const Function3D& f) {
    return 3;
}
std::size_t function_dimensions(const Function4D& f) {
    return 4;
}
std::size_t function_dimensions(const Function& f) {
    return std::visit([] (const auto& fp) { return function_dimensions(fp); },f); 
}

Function function_from_commandline(int argc, char **argv){
   	for (int i = 0; i < argc; ++i) {
		if (std::string(argv[i]) == "-function1d")
            return function1d(++i,argc,argv);
        else if (std::string(argv[i]) == "-function2d")
            return function2d(++i,argc,argv);
        else if (std::string(argv[i]) == "-function3d")
            return function3d(++i,argc,argv);
        else if (std::string(argv[i]) == "-function4d")
            return function4d(++i,argc,argv);
    }
    return Function1D(
        [] (double x) { return 0.5f; },
        [] (double a, double b) { return 0.5f*(b-a);});
}