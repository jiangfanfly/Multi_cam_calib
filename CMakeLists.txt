cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE Debug)
project(Muti_cam_calib)

IF(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)


find_package(OpenCV)
find_package(Eigen3 3.1.0 REQUIRED)
find_package(Ceres REQUIRED)

INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS}
${PROJECT_SOURCE_DIR}
${PROJECT_SOURCE_DIR}/include
)



add_library(${PROJECT_NAME} SHARED
src/exincalib.cc
)

MESSAGE("OpenCV_INCLUDE_DIRS: " ${OpenCV_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

link_directories(${OpenCV_LIBRARY_DIRS})

target_link_libraries(${PROJECT_NAME}
${OpenCV_LIBS}
${EIGEN3_LIBS}
${CERES_LIBRARIES}
)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJRCT_NAME}/bin)
add_executable(Multi_cam_calib_main "Multi_cam_calib_main.cpp.cpp")
target_link_libraries(Multi_cam_calib_main ${PROJECT_NAME})
