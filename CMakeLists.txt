project(plot)
cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

include(Paths) 
#  include(External)
include(Compiler)


###########################################################################################
# TARGETS
###########################################################################################
add_executable(test main/test.cc)

##########
# FOR DOCUMENTATION
##########
#add_executable(integrands main/doc/integrands.cc)
#add_executable(ranges main/doc/ranges.cc)
#add_executable(integrators main/doc/integrators.cc)
#add_executable(bins main/doc/bins.cc)
#add_executable(integrators-bins main/doc/integrators-bins.cc)

