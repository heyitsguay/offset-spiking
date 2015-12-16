cmake_minimum_required(VERSION 3.3)
project(offset_spiking)

set(SRCDIR, src)
set(BUILDDIR, build)
set(INCDIR, include)

set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -Wall -g")

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.55.0 COMPONENTS iostreams system filesystem REQUIRED)

find_package(OpenMP)

if(Boost_FOUND AND OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OpenMP_CXX_FLAGS}")

set(SOURCE_FILES src/main.cpp)
add_executable(offset_spiking ${SOURCE_FILES})