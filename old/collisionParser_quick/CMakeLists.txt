cmake_minimum_required(VERSION 3.10)
project(collisionParser_quick)

# Require 64-bit build
if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  message(FATAL_ERROR "This project requires a 64-bit compiler. Please use the x64 Native Tools Command Prompt for VS.")
endif()

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Optimization for MSVC
if(MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE "/O2")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif()

# ROOT
find_package(ROOT REQUIRED COMPONENTS Core Tree RIO Hist)
include(${ROOT_USE_FILE})

# ZLIB paths
set(ZLIB_BUILD_DIR "C:/Users/Sam/Documents/code/trimRunner/collisionParser_quick/zlib-1.3.1/build")
set(ZLIB_SOURCE_DIR "C:/Users/Sam/Documents/code/trimRunner/collisionParser_quick/zlib-1.3.1")
include_directories(${ZLIB_SOURCE_DIR} ${ZLIB_BUILD_DIR})
link_directories(${ZLIB_BUILD_DIR}/Release)

# Main executable
add_executable(collisionParser_quick collisionParser_quick.cpp gzstream.C)
target_include_directories(collisionParser_quick PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(collisionParser_quick PRIVATE ${ROOT_LIBRARIES} ${ZLIB_BUILD_DIR}/Release/zlibstatic.lib)
