# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_SOURCE_DIR = "/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code"

# Include any dependencies generated for this target.
include CMakeFiles/quad_mesh.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/quad_mesh.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/quad_mesh.dir/flags.make

CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o: CMakeFiles/quad_mesh.dir/flags.make
CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o: quad_mesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o -c "/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code/quad_mesh.cpp"

CMakeFiles/quad_mesh.dir/quad_mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quad_mesh.dir/quad_mesh.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code/quad_mesh.cpp" > CMakeFiles/quad_mesh.dir/quad_mesh.cpp.i

CMakeFiles/quad_mesh.dir/quad_mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quad_mesh.dir/quad_mesh.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code/quad_mesh.cpp" -o CMakeFiles/quad_mesh.dir/quad_mesh.cpp.s

CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o.requires:

.PHONY : CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o.requires

CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o.provides: CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o.requires
	$(MAKE) -f CMakeFiles/quad_mesh.dir/build.make CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o.provides.build
.PHONY : CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o.provides

CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o.provides.build: CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o


# Object files for target quad_mesh
quad_mesh_OBJECTS = \
"CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o"

# External object files for target quad_mesh
quad_mesh_EXTERNAL_OBJECTS =

quad_mesh: CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o
quad_mesh: CMakeFiles/quad_mesh.dir/build.make
quad_mesh: /usr/lib/x86_64-linux-gnu/libmpfr.so
quad_mesh: /usr/local/lib/libgmp.so
quad_mesh: /usr/local/lib/libCGAL.so.13.0.0
quad_mesh: /usr/local/lib/libboost_thread.so
quad_mesh: /usr/local/lib/libboost_system.so
quad_mesh: /usr/lib/x86_64-linux-gnu/libpthread.so
quad_mesh: /usr/local/lib/libCGAL.so.13.0.0
quad_mesh: /usr/local/lib/libboost_thread.so
quad_mesh: /usr/local/lib/libboost_system.so
quad_mesh: /usr/lib/x86_64-linux-gnu/libpthread.so
quad_mesh: /usr/lib/x86_64-linux-gnu/libglut.so
quad_mesh: /usr/lib/x86_64-linux-gnu/libXmu.so
quad_mesh: /usr/lib/x86_64-linux-gnu/libXi.so
quad_mesh: /usr/lib/x86_64-linux-gnu/libGLU.so
quad_mesh: /usr/lib/x86_64-linux-gnu/libGL.so
quad_mesh: CMakeFiles/quad_mesh.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable quad_mesh"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/quad_mesh.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/quad_mesh.dir/build: quad_mesh

.PHONY : CMakeFiles/quad_mesh.dir/build

CMakeFiles/quad_mesh.dir/requires: CMakeFiles/quad_mesh.dir/quad_mesh.cpp.o.requires

.PHONY : CMakeFiles/quad_mesh.dir/requires

CMakeFiles/quad_mesh.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/quad_mesh.dir/cmake_clean.cmake
.PHONY : CMakeFiles/quad_mesh.dir/clean

CMakeFiles/quad_mesh.dir/depend:
	cd "/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code" "/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code" "/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code" "/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code" "/home/amir/Documents/CS Masters/3dPrintingAlgorithms/Final Project/Code/CMakeFiles/quad_mesh.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/quad_mesh.dir/depend
