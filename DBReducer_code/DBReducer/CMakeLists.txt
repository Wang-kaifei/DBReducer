cmake_minimum_required(VERSION 3.15)
project(DBReducer)

set(CMAKE_CXX_STANDARD 14)
#SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

add_executable(DBReducer main.cpp Rebuild.cpp Param.cpp)
TARGET_LINK_LIBRARIES(DBReducer D:/DBReducer/DBReducer_code/utils/build/lib/libUtils.dll)