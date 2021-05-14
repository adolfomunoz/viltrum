###################################################################################
# COMPILER FLAGS
###################################################################################

# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE RelWithDebInfo)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_STANDARD 17)

option(BUILD_PROFILING "Build with profiling options" OFF)


# Select flags.
if ( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") )
    if (${BUILD_PROFILING})
        set(CMAKE_BUILD_TYPE RelWithDebInfo)
        set(CMAKE_CXX_FLAGS  "-Wall -pg")
    else()
        set(CMAKE_CXX_FLAGS  "-Wall")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -mtune=native")
    if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
       set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Dsrandom=srand -Drandom=rand -D_USE_MATH_DEFINES")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message("Using Visual Studio, are you sure?")
endif()

message(STATUS "Compiler  = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "System    = ${CMAKE_SYSTEM_NAME}")
message(STATUS "Prefix    = ${CMAKE_PREFIX_PATH}")
message(STATUS "Flags     = ${CMAKE_CXX_FLAGS}")
message(STATUS "Build     = ${CMAKE_BUILD_TYPE}")
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
   message(STATUS "R.Flags   = ${CMAKE_CXX_FLAGS_RELEASE}")
elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
   message(STATUS "D.Flags   = ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
   message(STATUS "D.Flags   = ${CMAKE_CXX_FLAGS_DEBUG}")
endif()


