# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ubuntu/program/github_rep/parallel_test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ubuntu/program/github_rep/parallel_test/build

# Include any dependencies generated for this target.
include CMakeFiles/my_executable.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/my_executable.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/my_executable.dir/flags.make

CMakeFiles/my_executable.dir/sources/test.cpp.o: CMakeFiles/my_executable.dir/flags.make
CMakeFiles/my_executable.dir/sources/test.cpp.o: ../sources/test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/program/github_rep/parallel_test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/my_executable.dir/sources/test.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_executable.dir/sources/test.cpp.o -c /home/ubuntu/program/github_rep/parallel_test/sources/test.cpp

CMakeFiles/my_executable.dir/sources/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_executable.dir/sources/test.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/program/github_rep/parallel_test/sources/test.cpp > CMakeFiles/my_executable.dir/sources/test.cpp.i

CMakeFiles/my_executable.dir/sources/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_executable.dir/sources/test.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/program/github_rep/parallel_test/sources/test.cpp -o CMakeFiles/my_executable.dir/sources/test.cpp.s

CMakeFiles/my_executable.dir/sources/test.cpp.o.requires:

.PHONY : CMakeFiles/my_executable.dir/sources/test.cpp.o.requires

CMakeFiles/my_executable.dir/sources/test.cpp.o.provides: CMakeFiles/my_executable.dir/sources/test.cpp.o.requires
	$(MAKE) -f CMakeFiles/my_executable.dir/build.make CMakeFiles/my_executable.dir/sources/test.cpp.o.provides.build
.PHONY : CMakeFiles/my_executable.dir/sources/test.cpp.o.provides

CMakeFiles/my_executable.dir/sources/test.cpp.o.provides.build: CMakeFiles/my_executable.dir/sources/test.cpp.o


CMakeFiles/my_executable.dir/sources/benchmark.cpp.o: CMakeFiles/my_executable.dir/flags.make
CMakeFiles/my_executable.dir/sources/benchmark.cpp.o: ../sources/benchmark.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/program/github_rep/parallel_test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/my_executable.dir/sources/benchmark.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_executable.dir/sources/benchmark.cpp.o -c /home/ubuntu/program/github_rep/parallel_test/sources/benchmark.cpp

CMakeFiles/my_executable.dir/sources/benchmark.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_executable.dir/sources/benchmark.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/program/github_rep/parallel_test/sources/benchmark.cpp > CMakeFiles/my_executable.dir/sources/benchmark.cpp.i

CMakeFiles/my_executable.dir/sources/benchmark.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_executable.dir/sources/benchmark.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/program/github_rep/parallel_test/sources/benchmark.cpp -o CMakeFiles/my_executable.dir/sources/benchmark.cpp.s

CMakeFiles/my_executable.dir/sources/benchmark.cpp.o.requires:

.PHONY : CMakeFiles/my_executable.dir/sources/benchmark.cpp.o.requires

CMakeFiles/my_executable.dir/sources/benchmark.cpp.o.provides: CMakeFiles/my_executable.dir/sources/benchmark.cpp.o.requires
	$(MAKE) -f CMakeFiles/my_executable.dir/build.make CMakeFiles/my_executable.dir/sources/benchmark.cpp.o.provides.build
.PHONY : CMakeFiles/my_executable.dir/sources/benchmark.cpp.o.provides

CMakeFiles/my_executable.dir/sources/benchmark.cpp.o.provides.build: CMakeFiles/my_executable.dir/sources/benchmark.cpp.o


# Object files for target my_executable
my_executable_OBJECTS = \
"CMakeFiles/my_executable.dir/sources/test.cpp.o" \
"CMakeFiles/my_executable.dir/sources/benchmark.cpp.o"

# External object files for target my_executable
my_executable_EXTERNAL_OBJECTS =

my_executable: CMakeFiles/my_executable.dir/sources/test.cpp.o
my_executable: CMakeFiles/my_executable.dir/sources/benchmark.cpp.o
my_executable: CMakeFiles/my_executable.dir/build.make
my_executable: libbenchmark.a
my_executable: /usr/lib/gcc/x86_64-linux-gnu/7/libgomp.so
my_executable: /usr/lib/x86_64-linux-gnu/libpthread.so
my_executable: CMakeFiles/my_executable.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/program/github_rep/parallel_test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable my_executable"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/my_executable.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/my_executable.dir/build: my_executable

.PHONY : CMakeFiles/my_executable.dir/build

CMakeFiles/my_executable.dir/requires: CMakeFiles/my_executable.dir/sources/test.cpp.o.requires
CMakeFiles/my_executable.dir/requires: CMakeFiles/my_executable.dir/sources/benchmark.cpp.o.requires

.PHONY : CMakeFiles/my_executable.dir/requires

CMakeFiles/my_executable.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/my_executable.dir/cmake_clean.cmake
.PHONY : CMakeFiles/my_executable.dir/clean

CMakeFiles/my_executable.dir/depend:
	cd /home/ubuntu/program/github_rep/parallel_test/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/program/github_rep/parallel_test /home/ubuntu/program/github_rep/parallel_test /home/ubuntu/program/github_rep/parallel_test/build /home/ubuntu/program/github_rep/parallel_test/build /home/ubuntu/program/github_rep/parallel_test/build/CMakeFiles/my_executable.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/my_executable.dir/depend

