# DBReducer 两个平台共用的编译配置（v0.4.0 起，与 pAnno 仓库的 CMakeLists.txt 同一套做法）
#
# 编译（推荐用 DBReducer_code/build.ps1 / build.sh）：
#   cmake -S <本目录> -B <编译目录>
#   cmake --build <编译目录> -j
# 产物：<编译目录>/bin/DBReducer（Windows 带 .exe）。utils 编成静态库直接链进程序，不再单独出 libUtils.dll。
# 两个平台都全静态（用户 11-10 定）：不依赖 MinGW 运行库 dll / 系统的 glibc、libstdc++、libgomp。

cmake_minimum_required(VERSION 3.16)
project(DBReducer LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(DBR_FULL_STATIC "完全静态链接（关掉时只静态链接 libstdc++/libgcc，开发调试用）" ON)
if(DBR_FULL_STATIC)
  set(DBR_LINK_STATIC -static)
else()
  set(DBR_LINK_STATIC -static-libstdc++ -static-libgcc)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(UTILS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../utils)

# utils（../utils/CMakeLists.txt 是早先单独编 libUtils.dll 用的，已不再使用）
add_library(dbr_utils STATIC ${UTILS_DIR}/IOUtil.cpp ${UTILS_DIR}/StringUtil.cpp)
target_compile_options(dbr_utils PRIVATE -O3 -Wall -std=gnu++14)

# 版本信息（每次编译刷新，写进状态文件的 tool_build 字段）
add_custom_target(dbr_build_info ALL
  COMMAND ${CMAKE_COMMAND} -DSRC=${CMAKE_CURRENT_SOURCE_DIR} -DOUT=${CMAKE_BINARY_DIR}/gen/dbr_build_info.h
          -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/build_info.cmake
  BYPRODUCTS ${CMAKE_BINARY_DIR}/gen/dbr_build_info.h
  COMMENT "版本信息")

# OpenMP 直接用 -fopenmp：find_package(OpenMP) 在 MinGW 上会链接 libgomp.dll.a，-static 也拦不住
add_executable(DBReducer main.cpp Rebuild.cpp Param.cpp Rebuild-spec.cpp Status.cpp ProcessGuard.cpp Precheck.cpp)
target_compile_options(DBReducer PRIVATE -fopenmp -O3 -Wall -std=gnu++20)
target_include_directories(DBReducer PRIVATE ${CMAKE_BINARY_DIR}/gen)
add_dependencies(DBReducer dbr_build_info)
target_link_options(DBReducer PRIVATE -fopenmp)
target_link_libraries(DBReducer PRIVATE dbr_utils ${DBR_LINK_STATIC})
