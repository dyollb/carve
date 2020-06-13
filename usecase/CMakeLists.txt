cmake_minimum_required(VERSION 2.8)

project(carve-preprocessing)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(VTK_DIR /Users/lloyd/Code/builds/vtk-8.2)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set(CARVE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include CACHE PATH "Carve library include dir")
set(CARVE_CONFIG_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include CACHE PATH "Carve library include dir")
set(CARVE_LIBRARIES "carve" CACHE STRING "Carve libraries")

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
