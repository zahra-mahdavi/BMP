cmake_minimum_required(VERSION 3.12)
project(bmp_refactor_v3 CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(bmp_headers INTERFACE)
target_include_directories(bmp_headers INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(bmp_demo src/demo_load.cpp)
target_link_libraries(bmp_demo PRIVATE bmp_headers)

add_executable(bmp_refactor src/main.cpp)
target_link_libraries(bmp_refactor PRIVATE bmp_headers)
