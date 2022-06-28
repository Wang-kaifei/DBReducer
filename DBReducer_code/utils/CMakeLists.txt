cmake_minimum_required(VERSION 3.20)
# 设置工程名称
PROJECT(PFINDLIB)
set(CMAKE_CXX_STANDARD 14)
#SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")
# 创建变量SRC 源码文件
SET(Utils IOUtil.cpp StringUtil.cpp)

SET(LIBRARY_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/lib)

ADD_LIBRARY(Utils SHARED  ${Utils})
#SET(LIBRARY_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/lib)
#SET(CMAKE_INSTALL_PREFIX "C:\\Users\\pFind\\Desktop\\libtest\\build")
# 安装动态链接库
#INSTALL(
 #   TARGETS  hello
    # 安装LIBRARY 动态链接库
  #  LIBRARY
#)

# 安装.h头文件
#INSTALL(FILES hello.h DESTINATION include/hello)