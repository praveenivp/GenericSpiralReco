cmake_minimum_required(VERSION 3.10)

project(MySpiralReco LANGUAGES CXX C)

find_package(ISMRMRD CONFIG REQUIRED)
find_package(Gadgetron REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_library(MySpiralReco SHARED
MyBucketToBufferGadget.cpp
RefRecoGadget.cpp
RemoveRefOSGadget.cpp
AddTrajInfoGadget.cpp
spiraltraj/vdspiral.cpp 
spiraltraj/spiraltraj_gadgetron.cpp
spiraltraj/nonCartesianTraj.cpp
SpiralGriddingGadget.cpp
)

target_link_libraries(MySpiralReco Gadgetron::gadgetron_core Gadgetron::gadgetron_toolbox_cpucore_math Gadgetron::gadgetron_toolbox_spiral Gadgetron::gadgetron_toolbox_cpunfft)