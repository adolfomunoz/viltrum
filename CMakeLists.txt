cmake_minimum_required(VERSION 3.10)
project(viltrum)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

include(Paths) 
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
target_link_libraries(test-region-sample PRIVATE TBB::tbb)
add_executable(test-region-sample-2d main/compilation-tests/region-sample-2d.cc)
target_link_libraries(test-region-sample-2d PRIVATE TBB::tbb)
add_executable(test-region-sample-3d main/compilation-tests/region-sample-3d.cc)
target_link_libraries(test-region-sample-3d PRIVATE TBB::tbb)
add_executable(montecarlo-rng main/compilation-tests/montecarlo-rng.cc)
target_link_libraries(montecarlo-rng PRIVATE TBB::tbb)


add_executable(crespo21 main/compilation-tests/crespo21.cc)
target_compile_options(crespo21 PRIVATE -pg)
target_link_options(crespo21 PRIVATE -pg)
target_link_libraries(crespo21 PRIVATE TBB::tbb)

##########
# FOR DOCUMENTATION
##########
add_executable(montecarlo-1d main/doc/montecarlo-1d.cc)
target_link_libraries(montecarlo-1d PRIVATE TBB::tbb)
add_executable(montecarlo-2d main/doc/montecarlo-2d.cc)
target_link_libraries(montecarlo-2d PRIVATE TBB::tbb)
add_executable(montecarlo-infd main/doc/montecarlo-infd.cc)
target_link_libraries(montecarlo-infd PRIVATE TBB::tbb)
add_executable(ranges main/doc/ranges.cc)
target_link_libraries(ranges PRIVATE TBB::tbb)
add_executable(ranges-infinite main/doc/ranges-infinite.cc)
target_link_libraries(ranges-infinite PRIVATE TBB::tbb)
add_executable(integrands main/doc/integrands.cc)
target_link_libraries(integrands PRIVATE TBB::tbb)
add_executable(integrands-infinite main/doc/integrands-infinite.cc)
target_link_libraries(integrands-infinite PRIVATE TBB::tbb)
add_executable(integrators main/doc/integrators.cc)
target_link_libraries(integrators PRIVATE TBB::tbb)
#add_executable(ranges main/doc/ranges.cc)
#add_executable(integrators main/doc/integrators.cc)
#add_executable(bins main/doc/bins.cc)
#add_executable(integrators-bins main/doc/integrators-bins.cc)

