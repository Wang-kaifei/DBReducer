# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\CMake\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\pFind\Desktop\pAnno_wkf\utils

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\pFind\Desktop\pAnno_wkf\utils\build

# Include any dependencies generated for this target.
include CMakeFiles/Utils.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Utils.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Utils.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Utils.dir/flags.make

CMakeFiles/Utils.dir/IOUtil.cpp.obj: CMakeFiles/Utils.dir/flags.make
CMakeFiles/Utils.dir/IOUtil.cpp.obj: ../IOUtil.cpp
CMakeFiles/Utils.dir/IOUtil.cpp.obj: CMakeFiles/Utils.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\pFind\Desktop\pAnno_wkf\utils\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Utils.dir/IOUtil.cpp.obj"
	C:\mingw-w64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Utils.dir/IOUtil.cpp.obj -MF CMakeFiles\Utils.dir\IOUtil.cpp.obj.d -o CMakeFiles\Utils.dir\IOUtil.cpp.obj -c C:\Users\pFind\Desktop\pAnno_wkf\utils\IOUtil.cpp

CMakeFiles/Utils.dir/IOUtil.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Utils.dir/IOUtil.cpp.i"
	C:\mingw-w64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\pFind\Desktop\pAnno_wkf\utils\IOUtil.cpp > CMakeFiles\Utils.dir\IOUtil.cpp.i

CMakeFiles/Utils.dir/IOUtil.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Utils.dir/IOUtil.cpp.s"
	C:\mingw-w64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\pFind\Desktop\pAnno_wkf\utils\IOUtil.cpp -o CMakeFiles\Utils.dir\IOUtil.cpp.s

CMakeFiles/Utils.dir/StringUtil.cpp.obj: CMakeFiles/Utils.dir/flags.make
CMakeFiles/Utils.dir/StringUtil.cpp.obj: ../StringUtil.cpp
CMakeFiles/Utils.dir/StringUtil.cpp.obj: CMakeFiles/Utils.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\pFind\Desktop\pAnno_wkf\utils\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Utils.dir/StringUtil.cpp.obj"
	C:\mingw-w64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Utils.dir/StringUtil.cpp.obj -MF CMakeFiles\Utils.dir\StringUtil.cpp.obj.d -o CMakeFiles\Utils.dir\StringUtil.cpp.obj -c C:\Users\pFind\Desktop\pAnno_wkf\utils\StringUtil.cpp

CMakeFiles/Utils.dir/StringUtil.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Utils.dir/StringUtil.cpp.i"
	C:\mingw-w64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\pFind\Desktop\pAnno_wkf\utils\StringUtil.cpp > CMakeFiles\Utils.dir\StringUtil.cpp.i

CMakeFiles/Utils.dir/StringUtil.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Utils.dir/StringUtil.cpp.s"
	C:\mingw-w64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\pFind\Desktop\pAnno_wkf\utils\StringUtil.cpp -o CMakeFiles\Utils.dir\StringUtil.cpp.s

# Object files for target Utils
Utils_OBJECTS = \
"CMakeFiles/Utils.dir/IOUtil.cpp.obj" \
"CMakeFiles/Utils.dir/StringUtil.cpp.obj"

# External object files for target Utils
Utils_EXTERNAL_OBJECTS =

lib/libUtils.dll: CMakeFiles/Utils.dir/IOUtil.cpp.obj
lib/libUtils.dll: CMakeFiles/Utils.dir/StringUtil.cpp.obj
lib/libUtils.dll: CMakeFiles/Utils.dir/build.make
lib/libUtils.dll: CMakeFiles/Utils.dir/linklibs.rsp
lib/libUtils.dll: CMakeFiles/Utils.dir/objects1.rsp
lib/libUtils.dll: CMakeFiles/Utils.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\pFind\Desktop\pAnno_wkf\utils\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX shared library lib\libUtils.dll"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Utils.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Utils.dir/build: lib/libUtils.dll
.PHONY : CMakeFiles/Utils.dir/build

CMakeFiles/Utils.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\Utils.dir\cmake_clean.cmake
.PHONY : CMakeFiles/Utils.dir/clean

CMakeFiles/Utils.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\pFind\Desktop\pAnno_wkf\utils C:\Users\pFind\Desktop\pAnno_wkf\utils C:\Users\pFind\Desktop\pAnno_wkf\utils\build C:\Users\pFind\Desktop\pAnno_wkf\utils\build C:\Users\pFind\Desktop\pAnno_wkf\utils\build\CMakeFiles\Utils.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Utils.dir/depend
