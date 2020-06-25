include_directories("${carve_SOURCE_DIR}/include")
include_directories("${carve_SOURCE_DIR}/external/GLOOP/include")
include_directories("${carve_SOURCE_DIR}/external/GLEW/include")
include_directories("${carve_SOURCE_DIR}/external/GLUI/include")

add_library(carve_fileformats
	read_ply.cpp
	read_ply.hpp
	write_ply.cpp
	write_ply.hpp
)
target_link_libraries(carve_fileformats	carve gloop_model)
set_target_properties(carve_fileformats PROPERTIES
    VERSION   "${carve_VERSION_MAJOR}.${carve_VERSION_MINOR}.${carve_VERSION_PATCH}"
    SOVERSION "${carve_VERSION_MAJOR}.${carve_VERSION_MINOR}"
    DEBUG_POSTFIX "_d")

add_library(carve_misc
	geometry.hpp
	geometry.cpp
)
target_link_libraries(carve_misc	carve)
set_target_properties(carve_misc PROPERTIES
    VERSION   "${carve_VERSION_MAJOR}.${carve_VERSION_MINOR}.${carve_VERSION_PATCH}"
    SOVERSION "${carve_VERSION_MAJOR}.${carve_VERSION_MINOR}"
    DEBUG_POSTFIX "_d")

if(CARVE_WITH_GUI)
  add_library(carve_ui             geom_draw.cpp scene.cpp)
endif(CARVE_WITH_GUI)
