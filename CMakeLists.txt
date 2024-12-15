cmake_minimum_required(VERSION 3.5.0)
project(vdpApp VERSION 0.1.0 LANGUAGES C CXX)

find_package(wxWidgets REQUIRED COMPONENTS net core base)
include(${wxWidgets_USE_FILE})
add_executable(vdpApp vdpApp.cpp)
target_link_libraries(vdpApp ${wxWidgets_LIBRARIES} mathplot)
