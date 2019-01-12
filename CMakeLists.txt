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
    message("Build without OpenCV")
endif(OpenCV_FOUND)

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
