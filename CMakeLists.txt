project(plot)
cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

include(Paths)
include(External)
include(Compiler)


###########################################################################################
# TARGETS
###########################################################################################
add_executable(test-multiarray-slice main/multiarray/test-slice.cc)
add_executable(test-multiarray-fold main/multiarray/test-fold.cc)
add_executable(test-multiarray-transform main/multiarray/test-transform.cc)
add_executable(test-multiarray-fill main/multiarray/test-fill.cc)
add_executable(test-multiarray-split main/multiarray/test-split.cc)
add_executable(test-region main/test-region.cc)
add_executable(test-quadrature main/test-quadrature.cc)
add_dependencies(test-quadrature ${function_1d_deps})
add_executable(test-quadrature-bins main/test-quadrature-bins.cc)
add_dependencies(test-quadrature-bins ${function_1d_deps})
#add_executable(test-quadrature-plot main/test-quadrature-plot.cc)
#add_dependencies(test-quadrature-plot svg-cpp-plot ${function_1d_deps})
add_executable(test-multidimensional-range main/test-multidimensional-range.cc)


#add_executable(convergence-1d main/convergence-1d.cc)
#add_dependencies(convergence-1d svg-cpp-plot ${function_1d_deps})
#target_link_libraries(convergence-1d ${svg_cpp_plot_libs} ${function_1d_libs})
#target_compile_definitions(convergence-1d PRIVATE ${svg_cpp_plot_defs} ${function_1d_defs})

#add_executable(convergence-2d main/convergence-2d.cc)
#add_dependencies(convergence-2d svg-cpp-plot ${function_2d_deps})
#target_link_libraries(convergence-2d ${svg_cpp_plot_libs} ${function_2d_libs})
#target_compile_definitions(convergence-2d PRIVATE ${svg_cpp_plot_defs} ${function_2d_defs})

#add_executable(convergence-project-2d main/convergence-project-2d.cc)
#add_dependencies(convergence-project-2d svg-cpp-plot ${function_2d_deps})
#target_link_libraries(convergence-project-2d ${svg_cpp_plot_libs} ${function_2d_libs})
#target_compile_definitions(convergence-project-2d PRIVATE ${svg_cpp_plot_defs} ${function_2d_defs})

#add_executable(convergence-controlvariates-1d main/convergence-controlvariates-1d.cc)
#add_dependencies(convergence-controlvariates-1d svg-cpp-plot ${function_1d_deps})
#target_link_libraries(convergence-controlvariates-1d ${svg_cpp_plot_libs} ${function_1d_libs})
#target_compile_definitions(convergence-controlvariates-1d PRIVATE ${svg_cpp_plot_defs} ${function_1d_defs})

add_executable(stepper-samples-count main/stepper-samples-count.cc)

#add_executable(transmittance_test tracking/transmittance_test.cpp)
#add_dependencies(transmittance_test svg-cpp-plot ${function_1d_deps})
#target_link_libraries(transmittance_test ${svg_cpp_plot_libs} ${function_1d_libs})
#target_compile_definitions(transmittance_test PRIVATE ${svg_cpp_plot_defs} ${function_1d_defs})

#add_executable(render4d main/render4d.cc)
#add_dependencies(render4d cimg cimg-additions mj2 eigen)
#target_link_libraries(render4d ${cimg_libs})
#target_compile_definitions(render4d PRIVATE ${cimg_defs})

#add_executable(render3d main/render3d.cc)
#add_dependencies(render3d cimg cimg-additions mj2 eigen)
#target_link_libraries(render3d ${cimg_libs})
#target_compile_definitions(render3d PRIVATE ${cimg_defs})

#add_executable(experiment-single-scattering main/experiment-single-scattering.cc)
#add_dependencies(experiment-single-scattering  cimg cimg-additions mj2 eigen)
#target_link_libraries(experiment-single-scattering ${cimg_libs})
#target_compile_definitions(experiment-single-scattering PRIVATE ${cimg_defs})

#add_executable(experiment-multiple-scattering main/experiment-multiple-scattering.cc)
#add_dependencies(experiment-multiple-scattering  cimg cimg-additions mj2 eigen)
#target_link_libraries(experiment-multiple-scattering ${cimg_libs})
#target_compile_definitions(experiment-multiple-scattering PRIVATE ${cimg_defs})

#add_executable(plot-project-2d main/plot-project-2d.cc)
#add_dependencies(plot-project-2d svg-cpp-plot ${function_2d_deps})
#target_link_libraries(plot-project-2d ${svg_cpp_plot_libs} ${function_2d_libs})
#target_compile_definitions(plot-project-2d PRIVATE ${svg_cpp_plot_defs} ${function_2d_defs})

add_executable(overview main/overview.cc)
add_dependencies(overview svg-cpp-plot ${function_1d_deps})
target_link_libraries(overview ${svg_cpp_plot_libs} ${function_1d_libs})
target_compile_definitions(overview PRIVATE ${svg_cpp_plot_defs} ${function_1d_defs})

#add_executable(plot-2d main/plot-2d.cc)
#add_dependencies(plot-2d svg-cpp-plot ${function_2d_deps})
#target_link_libraries(plot-2d ${svg_cpp_plot_libs} ${function_2d_libs})
#target_compile_definitions(plot-2d PRIVATE ${svg_cpp_plot_defs} ${function_2d_defs})

#add_executable(map-convergence-project-2d main/map-convergence-project-2d.cc)
#add_dependencies(map-convergence-project-2d svg-cpp-plot ${function_2d_deps})
#target_link_libraries(map-convergence-project-2d ${svg_cpp_plot_libs} ${function_2d_libs})
#target_compile_definitions(map-convergence-project-2d PRIVATE ${svg_cpp_plot_defs} ${function_2d_defs})

#add_executable(map-convergence-average main/map-convergence-average)
#add_dependencies(map-convergence-average svg-cpp-plot)
#target_link_libraries(map-convergence-average ${svg_cpp_plot_libs})
#target_compile_definitions(map-convergence-average PRIVATE ${svg_cpp_plot_defs})

#add_executable(plot main/plot.cc)
#add_dependencies(plot svg-cpp-plot ${function_2d_deps} ${function_1d_deps})
#target_link_libraries(plot ${svg_cpp_plot_libs} ${function_1d_libs} ${function_2d_libs})
#target_compile_definitions(plot PRIVATE ${svg_cpp_plot_defs} ${function_1d_defs} ${function_2d_defs})

#add_executable(map-convergence-project main/map-convergence-project.cc)
#add_dependencies(map-convergence-project svg-cpp-plot ${function_2d_deps} ${function_1d_deps})
#target_link_libraries(map-convergence-project ${svg_cpp_plot_libs} ${function_1d_libs} ${function_2d_libs})
#target_compile_definitions(map-convergence-project PRIVATE ${svg_cpp_plot_defs} ${function_1d_defs} ${function_2d_defs})

