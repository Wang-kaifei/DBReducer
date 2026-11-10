#!/bin/bash
# Linux 上编译 DBReducer（CMake，配置见 DBReducer/CMakeLists.txt）
#
# 用法:  ./build.sh <编译目录>
#        可用环境变量 CXX/CC 指定编译器（默认 g++-12/gcc-12，没有就用 g++/gcc）
# 产出:  <编译目录>/bin/DBReducer（全静态，CentOS 7 起都能跑）
set -e
OUT=${1:?用法: build.sh <编译目录>}
SRC=$(cd "$(dirname "$0")" && pwd)
if [ -z "$CXX" ]; then
    if command -v g++-12 >/dev/null; then CXX=g++-12; CC=gcc-12; else CXX=g++; CC=gcc; fi
fi
GEN=
command -v ninja >/dev/null && GEN="-G Ninja"
cmake -S "$SRC/DBReducer" -B "$OUT" $GEN -DCMAKE_C_COMPILER="${CC:-gcc}" -DCMAKE_CXX_COMPILER="$CXX"
cmake --build "$OUT" -j"$(nproc)"
ls -l "$OUT/bin"
