# 生成 dbr_build_info.h：#define DBR_BUILD "<日期> <git describe>"
# 用法：cmake -DSRC=<源码目录> -DOUT=<输出文件> -P build_info.cmake
# 内容没变就不重写，免得每次都重编 Status.cpp
execute_process(COMMAND git -C "${SRC}" describe --always --dirty
                OUTPUT_VARIABLE GIT_DESC OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT GIT_DESC)
  set(GIT_DESC "nogit")
endif()
string(TIMESTAMP TODAY "%Y-%m-%d")
set(CONTENT "#define DBR_BUILD \"${TODAY} ${GIT_DESC}\"\n")
if(EXISTS "${OUT}")
  file(READ "${OUT}" OLD)
else()
  set(OLD "")
endif()
if(NOT OLD STREQUAL CONTENT)
  file(WRITE "${OUT}" "${CONTENT}")
endif()
