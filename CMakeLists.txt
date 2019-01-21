# project config
cmake_minimum_required(VERSION 3.5)
project(masscollapse VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/externals/sanitizers-cmake/cmake" ${CMAKE_MODULE_PATH})
set(SRC src/Particle.cpp src/main.cpp)

find_package(Sanitizers)

# executable
add_executable(masscollapse ${SRC})
target_compile_features(masscollapse PRIVATE cxx_std_11)
set_target_properties(masscollapse PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_sanitizers(masscollapse)

# OpenCV
find_package(OpenCV QUIET)
if(OpenCV_FOUND)
    message("Build with OpenCV")
    target_compile_definitions(masscollapse PUBLIC USE_OPENCV)
    target_link_libraries(masscollapse "opencv_core" "opencv_videoio" "opencv_highgui")
else(OpenCV_FOUND)
    message("Build withOUT OpenCV")
endif(OpenCV_FOUND)

# OpenCL
include(FindOpenCL)
find_package(OpenCL REQUIRED)
include_directories(${OPENCL_INCLUDE_DIR})
if(OpenCL_FOUND)
	target_link_libraries(masscollapse ${OPENCL_LIBRARIES})
    target_compile_definitions(masscollapse PUBLIC USE_OPENCL)
    message("Build with OpenCL")
else(OpenCL_FOUND)
    message("Build withOUT OpenCL")
endif(OpenCL_FOUND)

# OpenCL Kernel file
configure_file(src/GravitationRK.cl GravitationRK.cl COPYONLY)

# OpenMP
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	find_package(OpenMP QUIET)
	if (OpenMP_FOUND)
	    message("Build with OpenMP")
	    target_compile_definitions(masscollapse PUBLIC USE_OPENMP)
	    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	else(OpenMP_FOUND)
	    message("Build withOUT OpenMP")
	endif(OpenMP_FOUND)
endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

# clang format
file(GLOB_RECURSE ALL_SOURCE_FILES src/*.cpp src/*.h)
add_custom_target(
    format
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND clang-format
    -style=file
    -i
    ${ALL_SOURCE_FILES}
)
