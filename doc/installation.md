# `viltrum` - Installation

`viltrum` is a header-only library, so installation is rather simple (it does not require any compilation process). You just need to download the library in a specific folder:

```
folder> git clone https://github.com/adolfomunoz/viltrum.git
```

Then make sure `folder` is within the included directories (`-Ifolder` parameter in g++, `include_directories("folder")` in CMake) and 
include it from C++.

```cpp
#include <viltrum/viltrum.h>
```

That would be enough to use all the features of the library.

There are other alternatives that might be more comfortable for you. If your project is already in git, `viltrum` can also be included as a submodule with git:

```
git submodule add https://github.com/adolfomunoz/viltrum.git
```

If you are using CMake as your compilation system, you can include `viltrum` as an external project:

```cmake
if (NOT EXTERNAL_INSTALL_LOCATION)
	set(EXTERNAL_INSTALL_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/external)
endif()
if (NOT IS_DIRECTORY ${EXTERNAL_INSTALL_LOCATION})
	file(MAKE_DIRECTORY ${EXTERNAL_INSTALL_LOCATION})
endif()

include(ExternalProject)
# External include directory
include_directories(${EXTERNAL_INSTALL_LOCATION})
ExternalProject_Add(viltrum
  GIT_REPOSITORY https://github.com/adolfomunoz/viltrum.git
  SOURCE_DIR ${EXTERNAL_INSTALL_LOCATION}/viltrum
  UPDATE_DISCONNECTED 1
  STEP_TARGETS update
  BUILD_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
)
```

and then for your `executable` that uses the library:

```cmake
add_executable(executable ${source_files})
add_dependencies(executable viltrum)
```

This will download `viltrum` right before compiling `executable`. If you need to update `viltrum` you can 
use the target `viltrum-update` as in:

```
make viltrum-update
```

or you can setup `UPDATE_DISCONNECTED 0` so it gets updated every time you compile (but this gets annoying quite soon).
