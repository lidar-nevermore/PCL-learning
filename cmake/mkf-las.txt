cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(lastest)

find_package(PCL 1.5 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS} $ENV{LASLIB_ROOT}/inc)
link_directories(${PCL_LIBRARY_DIRS} $ENV{LASLIB_ROOT}/lib)
add_definitions(${PCL_DEFINITIONS})

add_executable (lastest lastest.cpp)
target_link_libraries (lastest ${PCL_LIBRARIES} LASlibD.lib)
