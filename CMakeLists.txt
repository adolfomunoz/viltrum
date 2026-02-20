cmake_minimum_required(VERSION 3.10)
project(viltrum)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

include(Paths) 
#  include(External)
include(Compiler)


###########################################################################################
# TARGETS
###########################################################################################
find_package(TBB REQUIRED)

add_executable(test-array-parameter main/compilation-tests/array-parameter.cc)
target_link_libraries(test-array-parameter PRIVATE TBB::tbb)
add_executable(test-multiple-parameters main/compilation-tests/multiple-parameters.cc)
target_link_libraries(test-multiple-parameters PRIVATE TBB::tbb)
add_executable(test-multiple-parameters-2d main/compilation-tests/multiple-parameters-2d.cc)
target_link_libraries(test-multiple-parameters-2d PRIVATE TBB::tbb)
add_executable(test-multiple-parameters-sequence main/compilation-tests/multiple-parameters-sequence.cc)
target_link_libraries(test-multiple-parameters-sequence PRIVATE TBB::tbb)
add_executable(test-random-sequence main/compilation-tests/random-sequence.cc)
add_executable(test-concat main/compilation-tests/concat.cc)
add_executable(test-covariance main/compilation-tests/covariance.cc)
add_executable(test-inverse main/compilation-tests/inverse.cc)
add_executable(test-region-sample main/compilation-tests/region-sample.cc)
add_executable(test-region-sample-2d main/compilation-tests/region-sample-2d.cc)
add_executable(test-region-sample-3d main/compilation-tests/region-sample-3d.cc)

##########
# FOR DOCUMENTATION
##########
add_executable(montecarlo-1d main/doc/montecarlo-1d.cc)
add_executable(montecarlo-2d main/doc/montecarlo-2d.cc)
add_executable(montecarlo-infd main/doc/montecarlo-infd.cc)
add_executable(ranges main/doc/ranges.cc)
add_executable(ranges-infinite main/doc/ranges-infinite.cc)
add_executable(integrands main/doc/integrands.cc)
#add_executable(ranges main/doc/ranges.cc)
#add_executable(integrators main/doc/integrators.cc)
#add_executable(bins main/doc/bins.cc)
#add_executable(integrators-bins main/doc/integrators-bins.cc)

