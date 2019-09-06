cmake_minimum_required(VERSION 3.2.2)
project(Triangulation)
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -pedantic -Wextra -Wconversion -march=native -Wall -I. -ftest-coverage -fprofile-arcs")
find_package(OpenCV REQUIRED)

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

set (HEADERS "Poly.h" "PolyAbs.h" "PolyBase.h" "TriangulationBase.h" "LinearLS.h" "IterativeLS.h" "LinearEigen.h" "IterativeEigen.h")
set (SOURCES "Poly.cpp" "PolyAbs.cpp" "PolyBase.cpp" "TriangulationBase.cpp" "LinearLS.cpp" "IterativeLS.cpp" "LinearEigen.cpp" "IterativeEigen.cpp")

add_executable(Triangulation main.cpp ${SOURCES})

target_link_libraries(Triangulation ${OpenCV_LIBS})
target_link_libraries(Triangulation ${GTEST_LIBRARIES} pthread gcov)

add_test(AllTests Triangulation)

add_library(triangulation SHARED ${SOURCES})
set_target_properties(triangulation PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(triangulation PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS}")

install(TARGETS triangulation LIBRARY DESTINATION lib/Triangulation INCLUDES DESTINATION include/Triangulation)
install(FILES ${HEADERS} DESTINATION include/Triangulation)
configure_file(TriangulationConfig.cmake ${CMAKE_INSTALL_PREFIX}/lib/cmake/Triangulation/TriangulationConfig.cmake)

