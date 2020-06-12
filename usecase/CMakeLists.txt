cmake_minimum_required(VERSION 2.8)

project(carve-preprocessing)
SET(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin/)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set(CARVE_INCLUDE_DIR "" CACHE PATH "Carve library include dir")
set(CARVE_CONFIG_INCLUDE_DIR "" CACHE PATH "Carve library include dir")
set(CARVE_LIBRARIES "" CACHE STRING "Carve libraries")

include_directories(${CARVE_INCLUDE_DIR} ${CARVE_CONFIG_INCLUDE_DIR})

add_executable(prepro
    main.cpp
    VkiMultiDomainPreprocessing.cpp
)
target_link_libraries(prepro ${VTK_LIBRARIES} ${CARVE_LIBRARIES})

add_executable(boolean
    test_boolean.cpp
    Carve.cpp
)
target_link_libraries(boolean ${VTK_LIBRARIES} ${CARVE_LIBRARIES})
