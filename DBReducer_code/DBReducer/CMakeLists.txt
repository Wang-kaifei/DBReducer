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

# 2026-10：原先链接写死的 D:/DBReducer/DBReducer_code/utils/build/lib/libUtils.dll（他人目录），
# 改为把本仓库的 utils 作为子工程一起编译、按目标名链接。
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../utils ${CMAKE_CURRENT_BINARY_DIR}/utils)

# 版本信息写进 exe（状态文件的 tool_build 字段）
execute_process(COMMAND git describe --always --dirty
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_VARIABLE DBR_GIT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
string(TIMESTAMP DBR_DATE "%Y-%m-%d")
if(NOT DBR_GIT)
  set(DBR_GIT "nogit")
endif()

add_executable(DBReducer main.cpp Rebuild.cpp Param.cpp Rebuild-spec.cpp Status.cpp ProcessGuard.cpp Precheck.cpp)
target_compile_definitions(DBReducer PRIVATE DBR_BUILD="${DBR_DATE} ${DBR_GIT}")
TARGET_LINK_LIBRARIES(DBReducer Utils)
target_compile_options(DBReducer PUBLIC -fopenmp)