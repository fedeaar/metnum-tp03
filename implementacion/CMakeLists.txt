cmake_minimum_required(VERSION 3.16)
project(iter)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=gnu++14 -ggdb3 -g -O3 -DNDEBUG")

add_library(gtest STATIC dependencies/gtest-1.8.1/gtest-all.cc)
add_library(gtest_main STATIC dependencies/gtest-1.8.1/gtest_main.cc)
add_subdirectory(dependencies/eigen)

set(SOURCE_FILES
        src/pagerank.h
        src/impl/pagerank.cpp
        src/iter.h
        src/impl/iter.cpp
        src/grafo.h
        src/impl/grafo.cpp
        src/IO.h
        src/impl/IO.cpp)

add_executable(TESTS-iterativos-EG ${SOURCE_FILES} tests/iterativos-tests.cpp)
target_compile_definitions(TESTS-iterativos-EG PUBLIC METODO=0)
target_link_libraries(TESTS-iterativos-EG gtest gtest_main eigen)

add_executable(TESTS-iterativos-GS ${SOURCE_FILES} tests/iterativos-tests.cpp)
target_compile_definitions(TESTS-iterativos-GS PUBLIC METODO=1)
target_link_libraries(TESTS-iterativos-GS gtest gtest_main eigen)

add_executable(TESTS-iterativos-J ${SOURCE_FILES} tests/iterativos-tests.cpp)
target_compile_definitions(TESTS-iterativos-J PUBLIC METODO=2)
target_link_libraries(TESTS-iterativos-J gtest gtest_main eigen)

add_executable(tp3 ${SOURCE_FILES} src/run.h src/impl/run.cpp main.cpp)
target_link_libraries(tp3 eigen)