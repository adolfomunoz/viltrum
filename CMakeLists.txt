project(plot)
cmake_minimum_required(VERSION 3.10)

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
add_executable(test-random-sequence main/compilation-tests/random-sequence.cc)

##########
# FOR DOCUMENTATION
##########
#add_executable(integrands main/doc/integrands.cc)
#add_executable(ranges main/doc/ranges.cc)
#add_executable(integrators main/doc/integrators.cc)
#add_executable(bins main/doc/bins.cc)
#add_executable(integrators-bins main/doc/integrators-bins.cc)

