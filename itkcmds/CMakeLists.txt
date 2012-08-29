project(itkcmds)
cmake_minimum_required(VERSION 2.8)

find_package(ITK 4 REQUIRED COMPONENTS ITK-IO-Base)

add_executable(itkinfo itkinfo.cpp)
target_link_libraries(itkinfo ${ITK-IO-Base_Libraries})