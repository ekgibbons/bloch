set(LIBRARY_NAME "bloch_cardiac")

cmake_minimum_required(VERSION 2.8.12)
project(${LIBRARY_NAME})

# Set source directory
set(SOURCE_DIR "src")

# Set python directory
# set(PYTHON_DIR "${PROJECT_BINARY_DIR}/python")

# Tell CMake that headers are also in SOURCE_DIR
include_directories(${SOURCE_DIR})
file(GLOB SOURCES "src/*.cpp")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wconversion -O3")


# Pull in all of the Armadillo stuff
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Generate Python module
find_package(pybind11 REQUIRED)
pybind11_add_module(${LIBRARY_NAME} ${SOURCES})
# set_target_properties(${LIBRARY_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_link_libraries(${LIBRARY_NAME} PRIVATE ${ARMADILLO_LIBRARIES})



