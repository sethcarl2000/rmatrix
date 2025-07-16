#---------------------------------------------------------
#   I copied most of this CMakeLists.txt from: https://cliutils.gitlab.io/modern-cmake/examples/root-dict/README.html 
cmake_minimum_required(VERSION 3.11)

#---------------------------------------------------------
#   Define the project 
project(RMatrixLibrary VERSION 0.1
        LANGUAGES CXX)  

#---------------------------------------------------------
#   Which include directories to use. 
set(CMAKE_CXX_STANDARD
    20
    CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)


#---------------------------------------------------------
#   Find ROOT (necessary)
find_package(ROOT 6.20 REQUIRED COMPONENTS Core RIO)

#---------------------------------------------------------
#   Add Sources & Headers 

#---------------------------------------------------------
#   Generate the ROOT dictionary
root_generate_dictionary(G__RMatrix RMatrix.h LINKDEF LinkDef.h)

#---------------------------------------------------------
#   Add Libraries
add_library(RMatrix SHARED RMatrix.h RMatrix.cxx G__RMatrix.cxx)
target_include_directories(RMatrix PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(RMatrix PUBLIC ROOT::Core)
