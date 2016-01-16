cmake_minimum_required(VERSION 3.3)
project(offset-spiking)

include_directories(src include lib)

set(SRC
        src/KenyonCell.cpp
        src/Network.cpp
        src/CholinergicSynapse.cpp
        src/epsp_histogram.cpp
        src/main.cpp
        src/offset_test.cpp
        src/run_once.cpp
        src/input_test.cpp
        src/kc_threshold.cpp
        src/VoltageClamp.cpp
        src/plot_Vs.cpp
        src/clamp_test.cpp
)

set(INCLUDE
        include/KenyonCell.h
        include/Network.h
        include/CholinergicSynapse.h
        include/epsp_histogram.h
        include/offset_test.h
        include/run_once.h
        include/input_test.h
        include/kc_threshold.h
        include/time_constant.h
        include/VoltageClamp.h
        include/clamp_test.h
        include/plot_Vs.h
)


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

    add_executable(offset-spiking ${SRC} ${INCLUDE})
    target_link_libraries(offset-spiking ${Boost_LIBRARIES})
endif()