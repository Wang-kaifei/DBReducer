cmake_minimum_required(VERSION 3.15)
project(DBReducer)

set(CMAKE_CXX_STANDARD 20)
#SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")
FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(DBReducer main.cpp Rebuild.cpp Param.cpp Rebuild-spec.cpp)
TARGET_LINK_LIBRARIES(DBReducer D:/DBReducer/DBReducer_code/utils/build/lib/libUtils.dll)
target_compile_options(DBReducer PUBLIC -fopenmp)