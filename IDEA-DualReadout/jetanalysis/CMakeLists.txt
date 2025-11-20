cmake_minimum_required(VERSION 3.16)
project(jetResolution LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ROOT
find_package(ROOT REQUIRED COMPONENTS Core Hist Tree RIO Graf Gpad)
include(${ROOT_USE_FILE})

# PODIO + EDM4HEP
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(DD4hep REQUIRED)

# FASTJET from CVMFS
set(FASTJET_BASE "/cvmfs/sft.cern.ch/lcg/external/fastjet/3.0.3/x86_64-slc6-gcc48-opt")

include_directories(${FASTJET_BASE}/include)
link_directories(${FASTJET_BASE}/lib)

add_executable(jetResolution jetResolution.cpp)

target_link_libraries(jetResolution
    ROOT::Core ROOT::Hist ROOT::Tree ROOT::RIO ROOT::Graf ROOT::Gpad
    podio::podio
    podio::podioIO            # NECESSARIA PER makeReader()
    EDM4HEP::edm4hep
    DD4hep::DDCore
    DD4hep::DDRec
    DD4hep::DDG4
    fastjet
)

set_target_properties(jetResolution PROPERTIES
    BUILD_RPATH "${FASTJET_BASE}/lib"
)

