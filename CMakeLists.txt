cmake_minimum_required(VERSION 3.10)
project(refFFT LANGUAGES CXX C)


set(SPIRAL_HEADER_FILES RefFFT.h MyBucketToBufferGadget.h CFL_IO.h SpiralGriddingGadget.h)
set(SPIRAL_SRC_FILES  RefFFT.cpp MyBucketToBufferGadget.cpp SpiralGriddingGadget.cpp) 


find_package(ISMRMRD CONFIG REQUIRED)
find_package(Gadgetron REQUIRED)
find_package(cpr REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_library(refFFT SHARED
MyBucketToBufferGadget.cpp
RefRecoGadget.cpp
RemoveRefOSGadget.cpp
AddTrajInfoGadget.cpp
spiraltraj/vdspiral.cpp 
spiraltraj/spiraltraj_gadgetron.cpp
spiraltraj/nonCartesianTraj.cpp
SpiralGriddingGadget.cpp
)
#target_include_directories(refFFT PRIVATE /home/gadgetron/core)
target_link_libraries(refFFT Gadgetron::gadgetron_core Gadgetron::gadgetron_toolbox_cpucore_math Gadgetron::gadgetron_toolbox_spiral Gadgetron::gadgetron_toolbox_cpunfft)