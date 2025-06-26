add_library(quickhull OBJECT quickhull/QuickHull.cpp)
target_include_directories(quickhull PUBLIC
    $<INSTALL_INTERFACE:quickhull>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/quickhull>)

target_compile_features(quickhull PUBLIC cxx_std_17)
install(TARGETS quickhull EXPORT manifoldTargets)

set(VORO_BUILD_SHARED_LIBS OFF CACHE BOOL "")
set(VORO_BUILD_EXAMPLES    OFF CACHE BOOL "")
set(VORO_BUILD_CMD_LINE    OFF CACHE BOOL "")
set(VORO_ENABLE_DOXYGEN    OFF CACHE BOOL "")
add_subdirectory(voro)

# Configure CDT as a library
# CDT is set up as an executable, so we need to create a library from its sources
file(GLOB CDT_SOURCES 
  CDT/src/giftWrap.cpp
  CDT/src/PLC.cpp
  CDT/src/delaunay.cpp
  CDT/src/main.cpp
)
add_library(CDT STATIC ${CDT_SOURCES})
target_compile_features(CDT PUBLIC cxx_std_20)
# Add CDT compile options
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  target_compile_options(CDT PUBLIC /fp:strict /Oi "/D _CRT_SECURE_NO_WARNINGS" "/D USE_MAROTS_METHOD" "/D LGPL=ON")
else()
  target_compile_options(CDT PUBLIC -O2 -frounding-math "-DUSE_MAROTS_METHOD" "-DLGPL=ON")
endif()

target_include_directories(CDT PUBLIC
  $<INSTALL_INTERFACE:CDT>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CDT/src>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CDT/include>
)

install(TARGETS CDT EXPORT manifoldTargets)
