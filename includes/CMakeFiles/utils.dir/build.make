# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/alexthem/Desktop/project-hw3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/alexthem/Desktop/project-hw3

# Include any dependencies generated for this target.
include includes/CMakeFiles/utils.dir/depend.make

# Include the progress variables for this target.
include includes/CMakeFiles/utils.dir/progress.make

# Include the compile flags for this target's objects.
include includes/CMakeFiles/utils.dir/flags.make

includes/CMakeFiles/utils.dir/convex_hull_alg.cpp.o: includes/CMakeFiles/utils.dir/flags.make
includes/CMakeFiles/utils.dir/convex_hull_alg.cpp.o: includes/convex_hull_alg.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alexthem/Desktop/project-hw3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object includes/CMakeFiles/utils.dir/convex_hull_alg.cpp.o"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/convex_hull_alg.cpp.o -c /home/alexthem/Desktop/project-hw3/includes/convex_hull_alg.cpp

includes/CMakeFiles/utils.dir/convex_hull_alg.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/convex_hull_alg.cpp.i"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alexthem/Desktop/project-hw3/includes/convex_hull_alg.cpp > CMakeFiles/utils.dir/convex_hull_alg.cpp.i

includes/CMakeFiles/utils.dir/convex_hull_alg.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/convex_hull_alg.cpp.s"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alexthem/Desktop/project-hw3/includes/convex_hull_alg.cpp -o CMakeFiles/utils.dir/convex_hull_alg.cpp.s

includes/CMakeFiles/utils.dir/functions.cpp.o: includes/CMakeFiles/utils.dir/flags.make
includes/CMakeFiles/utils.dir/functions.cpp.o: includes/functions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alexthem/Desktop/project-hw3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object includes/CMakeFiles/utils.dir/functions.cpp.o"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/functions.cpp.o -c /home/alexthem/Desktop/project-hw3/includes/functions.cpp

includes/CMakeFiles/utils.dir/functions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/functions.cpp.i"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alexthem/Desktop/project-hw3/includes/functions.cpp > CMakeFiles/utils.dir/functions.cpp.i

includes/CMakeFiles/utils.dir/functions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/functions.cpp.s"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alexthem/Desktop/project-hw3/includes/functions.cpp -o CMakeFiles/utils.dir/functions.cpp.s

includes/CMakeFiles/utils.dir/incremental_alg.cpp.o: includes/CMakeFiles/utils.dir/flags.make
includes/CMakeFiles/utils.dir/incremental_alg.cpp.o: includes/incremental_alg.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alexthem/Desktop/project-hw3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object includes/CMakeFiles/utils.dir/incremental_alg.cpp.o"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/incremental_alg.cpp.o -c /home/alexthem/Desktop/project-hw3/includes/incremental_alg.cpp

includes/CMakeFiles/utils.dir/incremental_alg.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/incremental_alg.cpp.i"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alexthem/Desktop/project-hw3/includes/incremental_alg.cpp > CMakeFiles/utils.dir/incremental_alg.cpp.i

includes/CMakeFiles/utils.dir/incremental_alg.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/incremental_alg.cpp.s"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alexthem/Desktop/project-hw3/includes/incremental_alg.cpp -o CMakeFiles/utils.dir/incremental_alg.cpp.s

includes/CMakeFiles/utils.dir/incremental_alg_subdivisial.cpp.o: includes/CMakeFiles/utils.dir/flags.make
includes/CMakeFiles/utils.dir/incremental_alg_subdivisial.cpp.o: includes/incremental_alg_subdivisial.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alexthem/Desktop/project-hw3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object includes/CMakeFiles/utils.dir/incremental_alg_subdivisial.cpp.o"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/incremental_alg_subdivisial.cpp.o -c /home/alexthem/Desktop/project-hw3/includes/incremental_alg_subdivisial.cpp

includes/CMakeFiles/utils.dir/incremental_alg_subdivisial.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/incremental_alg_subdivisial.cpp.i"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alexthem/Desktop/project-hw3/includes/incremental_alg_subdivisial.cpp > CMakeFiles/utils.dir/incremental_alg_subdivisial.cpp.i

includes/CMakeFiles/utils.dir/incremental_alg_subdivisial.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/incremental_alg_subdivisial.cpp.s"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alexthem/Desktop/project-hw3/includes/incremental_alg_subdivisial.cpp -o CMakeFiles/utils.dir/incremental_alg_subdivisial.cpp.s

includes/CMakeFiles/utils.dir/localSearch.cpp.o: includes/CMakeFiles/utils.dir/flags.make
includes/CMakeFiles/utils.dir/localSearch.cpp.o: includes/localSearch.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alexthem/Desktop/project-hw3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object includes/CMakeFiles/utils.dir/localSearch.cpp.o"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/localSearch.cpp.o -c /home/alexthem/Desktop/project-hw3/includes/localSearch.cpp

includes/CMakeFiles/utils.dir/localSearch.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/localSearch.cpp.i"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alexthem/Desktop/project-hw3/includes/localSearch.cpp > CMakeFiles/utils.dir/localSearch.cpp.i

includes/CMakeFiles/utils.dir/localSearch.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/localSearch.cpp.s"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alexthem/Desktop/project-hw3/includes/localSearch.cpp -o CMakeFiles/utils.dir/localSearch.cpp.s

includes/CMakeFiles/utils.dir/simulatedAnnealing.cpp.o: includes/CMakeFiles/utils.dir/flags.make
includes/CMakeFiles/utils.dir/simulatedAnnealing.cpp.o: includes/simulatedAnnealing.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alexthem/Desktop/project-hw3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object includes/CMakeFiles/utils.dir/simulatedAnnealing.cpp.o"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/simulatedAnnealing.cpp.o -c /home/alexthem/Desktop/project-hw3/includes/simulatedAnnealing.cpp

includes/CMakeFiles/utils.dir/simulatedAnnealing.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/simulatedAnnealing.cpp.i"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alexthem/Desktop/project-hw3/includes/simulatedAnnealing.cpp > CMakeFiles/utils.dir/simulatedAnnealing.cpp.i

includes/CMakeFiles/utils.dir/simulatedAnnealing.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/simulatedAnnealing.cpp.s"
	cd /home/alexthem/Desktop/project-hw3/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alexthem/Desktop/project-hw3/includes/simulatedAnnealing.cpp -o CMakeFiles/utils.dir/simulatedAnnealing.cpp.s

# Object files for target utils
utils_OBJECTS = \
"CMakeFiles/utils.dir/convex_hull_alg.cpp.o" \
"CMakeFiles/utils.dir/functions.cpp.o" \
"CMakeFiles/utils.dir/incremental_alg.cpp.o" \
"CMakeFiles/utils.dir/incremental_alg_subdivisial.cpp.o" \
"CMakeFiles/utils.dir/localSearch.cpp.o" \
"CMakeFiles/utils.dir/simulatedAnnealing.cpp.o"

# External object files for target utils
utils_EXTERNAL_OBJECTS =

includes/libutils.a: includes/CMakeFiles/utils.dir/convex_hull_alg.cpp.o
includes/libutils.a: includes/CMakeFiles/utils.dir/functions.cpp.o
includes/libutils.a: includes/CMakeFiles/utils.dir/incremental_alg.cpp.o
includes/libutils.a: includes/CMakeFiles/utils.dir/incremental_alg_subdivisial.cpp.o
includes/libutils.a: includes/CMakeFiles/utils.dir/localSearch.cpp.o
includes/libutils.a: includes/CMakeFiles/utils.dir/simulatedAnnealing.cpp.o
includes/libutils.a: includes/CMakeFiles/utils.dir/build.make
includes/libutils.a: includes/CMakeFiles/utils.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/alexthem/Desktop/project-hw3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX static library libutils.a"
	cd /home/alexthem/Desktop/project-hw3/includes && $(CMAKE_COMMAND) -P CMakeFiles/utils.dir/cmake_clean_target.cmake
	cd /home/alexthem/Desktop/project-hw3/includes && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/utils.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
includes/CMakeFiles/utils.dir/build: includes/libutils.a

.PHONY : includes/CMakeFiles/utils.dir/build

includes/CMakeFiles/utils.dir/clean:
	cd /home/alexthem/Desktop/project-hw3/includes && $(CMAKE_COMMAND) -P CMakeFiles/utils.dir/cmake_clean.cmake
.PHONY : includes/CMakeFiles/utils.dir/clean

includes/CMakeFiles/utils.dir/depend:
	cd /home/alexthem/Desktop/project-hw3 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/alexthem/Desktop/project-hw3 /home/alexthem/Desktop/project-hw3/includes /home/alexthem/Desktop/project-hw3 /home/alexthem/Desktop/project-hw3/includes /home/alexthem/Desktop/project-hw3/includes/CMakeFiles/utils.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : includes/CMakeFiles/utils.dir/depend

