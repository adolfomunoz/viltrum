##################################################################################
# LIBPNG
##################################################################################
find_package(PNG)
if(PNG_FOUND)
	include_directories(${PNG_INCLUDE_DIR})
	list(APPEND cimg_libs ${PNG_LIBRARY})
	list(APPEND cimg_defs "cimg_use_png")
	list(APPEND svg_cpp_plot_libs ${PNG_LIBRARY})
	list(APPEND svg_cpp_plot_defs "USE_PNG")
endif(PNG_FOUND)

find_package(JPEG)
if(JPEG_FOUND)
	include_directories(${JPEG_INCLUDE_DIR})
	list(APPEND cimg_libs ${JPEG_LIBRARY})
	list(APPEND cimg_defs "cimg_use_jpeg")
endif(JPEG_FOUND)

find_package(TIFF)
if(TIFF_FOUND)
	include_directories(${TIFF_INCLUDE_DIR})
	list(APPEND cimg_libs ${TIFF_LIBRARY})
	list(APPEND image_defs "cimg_use_tiff")
endif(TIFF_FOUND)

if(FFTW3_FOUND)
	include_directories(${FFTW3_INCLUDE_DIR})
	list(APPEND cimg_libs ${FFTW3_LIBRARY})
	list(APPEND cimg_defs "cimg_use_fftw3")
endif(FFTW3_FOUND)

#find_package(ImageMagick COMPONENTS Magick++)
#if (ImageMagick_FOUND)
#	include_directories(${ImageMagick_INCLUDE_DIRS})
#	list(APPEND cimg_libs ${ImageMagick_LIBRARIES})
#	list(APPEND cimg_defs "cimg_use_magick")
#	list(APPEND cimg_defs "MAGICKCORE_QUANTUM_DEPTH=16")
#	list(APPEND cimg_defs "MAGICKCORE_HDRI_ENABLE=0")
#endif(ImageMagick_FOUND)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    # Windows specific code
    list(APPEND cimgdisplay_libs gdi32)
    list(APPEND cimgdisplay_defs cimg_display=2)
    list(APPEND cimg_defs cimg_OS=2)
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    list(APPEND cimg_libs pthread)
    list(APPEND cimgdisplay_libs X11)
    list(APPEND cimgdisplay_defs cimg_display=1)
    list(APPEND cimg_defs cimg_OS=1)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")


######################################################################
# EXTERNAL LIBRARIES
######################################################################
if (NOT EXTERNAL_INSTALL_LOCATION)
	set(EXTERNAL_INSTALL_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/external)
endif()
if (NOT IS_DIRECTORY ${EXTERNAL_INSTALL_LOCATION})
	file(MAKE_DIRECTORY ${EXTERNAL_INSTALL_LOCATION})
endif()


include(ExternalProject)
# External include directory
include_directories(${EXTERNAL_INSTALL_LOCATION})
add_custom_target(update)

ExternalProject_Add(svg-cpp-plot
  GIT_REPOSITORY https://github.com/adolfomunoz/svg-cpp-plot.git
  SOURCE_DIR ${EXTERNAL_INSTALL_LOCATION}/svg-cpp-plot
  UPDATE_DISCONNECTED 1
  STEP_TARGETS update
  BUILD_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
)
add_dependencies(update svg-cpp-plot-update)

ExternalProject_Add(perlin-noise
  GIT_REPOSITORY https://github.com/Reputeless/PerlinNoise
  SOURCE_DIR ${EXTERNAL_INSTALL_LOCATION}/PerlinNoise
  UPDATE_DISCONNECTED 1
  STEP_TARGETS update
  BUILD_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
)
add_dependencies(update perlin-noise-update)

ExternalProject_Add(cimg
  GIT_REPOSITORY https://github.com/dtschump/CImg.git 
  SOURCE_DIR ${EXTERNAL_INSTALL_LOCATION}/CImg
  UPDATE_DISCONNECTED 1
  STEP_TARGETS update
  BUILD_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
)
add_dependencies(update cimg-update)

ExternalProject_Add(cimg-additions
  GIT_REPOSITORY https://github.com/adolfomunoz/cimg_additions.git
  SOURCE_DIR ${EXTERNAL_INSTALL_LOCATION}/cimg_additions
  UPDATE_DISCONNECTED 1
  STEP_TARGETS update
  BUILD_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
)
add_dependencies(update cimg-additions-update)

ExternalProject_Add(mj2
  GIT_REPOSITORY https://github.com/adolfomunoz/mj2.git
  SOURCE_DIR ${EXTERNAL_INSTALL_LOCATION}/mj2
  UPDATE_DISCONNECTED 1
  STEP_TARGETS update
  BUILD_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
)
add_dependencies(update mj2-update)

ExternalProject_Add(eigen
  GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror
  SOURCE_DIR ${EXTERNAL_INSTALL_LOCATION}/eigen
  UPDATE_DISCONNECTED 1
  STEP_TARGETS update
  BUILD_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
)
add_dependencies(update eigen-update)



include_directories(${EXTERNAL_INSTALL_LOCATION}/CImg)
include_directories(${EXTERNAL_INSTALL_LOCATION}/cimg_additions)
include_directories(${EXTERNAL_INSTALL_LOCATION}/eigen)

set(function_2d_defs ${cimg_defs})
set(function_2d_libs ${cimg_libs})
set(function_2d_deps perlin-noise cimg cimg-additions)
set(function_1d_defs "")
set(function_1d_libs "")
set(function_1d_deps perlin-noise)






