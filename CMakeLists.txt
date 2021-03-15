cmake_minimum_required(VERSION 3.0)
project(SZ3)
#project(First VERSION 3.0.0.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-std=c++11 -g -O3")
include(GNUInstallDirs)

find_library(ZSTD_LIBS zstd HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ/install/lib")
set(ZSTD_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ/install/include")

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(
        ${PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_FULL_INCLUDEDIR}>
        ${ZSTD_INCLUDES})


install(DIRECTORY ${${PROJECT_NAME}_SOURCE_DIR}/include/  DESTINATION ${CMAKE_INSTALL_FULL_INCLUDEDIR})

add_subdirectory(test)

