cmake_minimum_required(VERSION 3.10)
project(collisionParser)

# Require 64-bit build
if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  message(FATAL_ERROR "This project requires a 64-bit compiler. Please use the x64 Native Tools Command Prompt for VS.")
endif()

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Core Tree RIO Hist)

# Include ROOT's headers
include(${ROOT_USE_FILE})

# Add your executable
add_executable(collisionParser collisionParser.cpp)

# Link against ROOT libraries
target_link_libraries(collisionParser ${ROOT_LIBRARIES})
