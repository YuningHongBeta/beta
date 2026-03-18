cmake_minimum_required(VERSION 3.16...3.27)
project(beta)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
add_executable(beta ${sources})
target_link_libraries(beta ${Geant4_LIBRARIES})

file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})