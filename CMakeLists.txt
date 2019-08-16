cmake_minimum_required(VERSION 3.2.2)
project(Triangulation)
set(CMAKE_CXX_STANDARD 11)
find_package(OpenCV REQUIRED)
find_package(GSL REQUIRED)

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(Triangulation main.cpp Poly.cpp PolyAbs.cpp PolyBase.cpp)

target_link_libraries(Triangulation ${OpenCV_LIBS})
target_link_libraries(Triangulation ${GSL_LIBRARY})
target_link_libraries(Triangulation ${GTEST_LIBRARIES} pthread gcov)

add_test(AllTests Triangulation)
