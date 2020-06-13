cmake_minimum_required(VERSION 3.0)

project(carve)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(carve_VERSION_MAJOR 2)
set(carve_VERSION_MINOR 1)
set(carve_VERSION_PATCH 0)

set(CARVE_VERSION ${carve_VERSION_MAJOR}.${carve_VERSION_MINOR}.${carve_VERSION_PATCH})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

option(BUILD_COVERAGE                    "Compile with gcov"                                 OFF)
option(BUILD_SHARED_LIBS                 "Compile libcarve as shared library/dll"            ON)
option(CARVE_WITH_GUI                    "Compile gui code"                                  OFF)
option(CARVE_DEBUG                       "Compile in debug code"                             OFF)
option(CARVE_DEBUG_WRITE_PLY_DATA        "Write geometry output during debug"                OFF)
option(CARVE_USE_EXACT_PREDICATES        "Use Shewchuk's exact predicates, where possible"   ON)
option(CARVE_INTERSECT_GLU_TRIANGULATOR  "Include support for GLU triangulator in intersect" OFF)
option(CARVE_GTEST_TESTS                 "Complie gtest, and dependent tests"                OFF)

if(MSVC)
  add_compile_options(/W3 /D_SCL_SECURE_NO_WARNINGS)# /WX)
else()
  add_compile_options(-Wall)# -Wextra -pedantic -Werror)
endif()

if (MSVC AND NOT BUILD_SHARED_LIBS)
  # For MSVC, CMake sets certain flags to defaults we want to override.
  # This replacement code is taken from sample in the CMake Wiki at
  # http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace.
  foreach (flag_var
           CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
           CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    # In hermetic build environments, tests may not have access to MS runtime
    # DLLs, so this replaces /MD (CRT libraries in DLLs) with /MT (static CRT
    # libraries).
    string(REPLACE "/MD" "-MT" ${flag_var} "${${flag_var}}")
    # We prefer more strict warning checking for building Google Test.
    # Replaces /W3 with /W4 in defaults.
    string(REPLACE "/W3" "-W4" ${flag_var} "${${flag_var}}")
  endforeach()
endif()

SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

if (BUILD_COVERAGE)
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_BUILD_TYPE DEBUG)
  IF (CMAKE_COMPILER_IS_GNUCC)
    SET(CMAKE_C_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
    SET(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -lgcov")
  ENDIF()
endif(BUILD_COVERAGE)

if(WIN32)
  if(BUILD_SHARED_LIBS)
    message(WARNING "BUILD_SHARED_LIBS=ON: not all symbols are exported yet")
  else(BUILD_SHARED_LIBS)
    add_definitions(-DCARVE_STATIC)
  endif(BUILD_SHARED_LIBS)
  add_definitions(-D_USE_MATH_DEFINES)
  add_definitions(-DNOMINMAX)
endif(WIN32)

find_package(OpenMP)

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(CARVE_WITH_GUI)
  find_package(OpenGL)
  find_package(GLUT)
  find_package(GLEW)

  if(NOT OPENGL_FOUND)
    message(WARNING "Unable to locate OpenGL")
    set(CARVE_WITH_GUI OFF)

  elseif(NOT GLUT_FOUND)
    message(WARNING "Unable to locate GLUT")
    set(CARVE_WITH_GUI OFF)

  else(OPENGL_FOUND AND GLUT_FOUND)
    message(STATUS "Found OpenGL and GLUT")
    include_directories(${OPENGL_INCLUDE_DIR})
    include_directories(${GLUT_INCLUDE_DIR})

    if(GLEW_FOUND)
      include_directories(${GLEW_INCLUDE_DIR})
    else(GLEW_FOUND)
      if(WIN32)
        add_definitions(-DGLEW_STATIC)
      endif(WIN32)
      add_subdirectory(external/GLEW) 
    endif(GLEW_FOUND)

    if(WIN32)
      add_definitions(-DGLUI_NO_LIB_PRAGMA)
      add_definitions(-DGLUI_USE_STATIC_LIB)
    endif(WIN32)
    add_subdirectory(external/GLUI) 

  endif(NOT OPENGL_FOUND)

endif(CARVE_WITH_GUI)

add_subdirectory(external/GLOOP) 

if (CARVE_GTEST_TESTS)
  add_subdirectory(external/googletest)
endif(CARVE_GTEST_TESTS)

configure_file (
  "${carve_SOURCE_DIR}/include/carve/cmake-config.h.in"
  "${carve_BINARY_DIR}/include/carve/config.h"
)
add_definitions(-DCMAKE_BUILD)

include_directories(${carve_BINARY_DIR}/include)

add_subdirectory(lib) 
add_subdirectory(include) 
add_subdirectory(common) 
add_subdirectory(src) 
add_subdirectory(examples) 
add_subdirectory(tests)

add_subdirectory(usecase)

include(CTest)

if(BUILD_TESTING)
  ENABLE_TESTING()
endif(BUILD_TESTING)
