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
CMAKE_SOURCE_DIR = /home/divij/Documents/ParMooN_Codes/ParMooN_CMG

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD

# Include any dependencies generated for this target.
include src/AMG/CMakeFiles/amg.dir/depend.make

# Include the progress variables for this target.
include src/AMG/CMakeFiles/amg.dir/progress.make

# Include the compile flags for this target's objects.
include src/AMG/CMakeFiles/amg.dir/flags.make

src/AMG/CMakeFiles/amg.dir/amg_1d_prec.c.o: src/AMG/CMakeFiles/amg.dir/flags.make
src/AMG/CMakeFiles/amg.dir/amg_1d_prec.c.o: ../src/AMG/amg_1d_prec.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/AMG/CMakeFiles/amg.dir/amg_1d_prec.c.o"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/amg.dir/amg_1d_prec.c.o   -c /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_1d_prec.c

src/AMG/CMakeFiles/amg.dir/amg_1d_prec.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/amg.dir/amg_1d_prec.c.i"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_1d_prec.c > CMakeFiles/amg.dir/amg_1d_prec.c.i

src/AMG/CMakeFiles/amg.dir/amg_1d_prec.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/amg.dir/amg_1d_prec.c.s"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_1d_prec.c -o CMakeFiles/amg.dir/amg_1d_prec.c.s

src/AMG/CMakeFiles/amg.dir/amg_2d_prec.c.o: src/AMG/CMakeFiles/amg.dir/flags.make
src/AMG/CMakeFiles/amg.dir/amg_2d_prec.c.o: ../src/AMG/amg_2d_prec.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/AMG/CMakeFiles/amg.dir/amg_2d_prec.c.o"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/amg.dir/amg_2d_prec.c.o   -c /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_2d_prec.c

src/AMG/CMakeFiles/amg.dir/amg_2d_prec.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/amg.dir/amg_2d_prec.c.i"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_2d_prec.c > CMakeFiles/amg.dir/amg_2d_prec.c.i

src/AMG/CMakeFiles/amg.dir/amg_2d_prec.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/amg.dir/amg_2d_prec.c.s"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_2d_prec.c -o CMakeFiles/amg.dir/amg_2d_prec.c.s

src/AMG/CMakeFiles/amg.dir/amg_blas.c.o: src/AMG/CMakeFiles/amg.dir/flags.make
src/AMG/CMakeFiles/amg.dir/amg_blas.c.o: ../src/AMG/amg_blas.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object src/AMG/CMakeFiles/amg.dir/amg_blas.c.o"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/amg.dir/amg_blas.c.o   -c /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_blas.c

src/AMG/CMakeFiles/amg.dir/amg_blas.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/amg.dir/amg_blas.c.i"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_blas.c > CMakeFiles/amg.dir/amg_blas.c.i

src/AMG/CMakeFiles/amg.dir/amg_blas.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/amg.dir/amg_blas.c.s"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_blas.c -o CMakeFiles/amg.dir/amg_blas.c.s

src/AMG/CMakeFiles/amg.dir/amg_blas_nse2d.c.o: src/AMG/CMakeFiles/amg.dir/flags.make
src/AMG/CMakeFiles/amg.dir/amg_blas_nse2d.c.o: ../src/AMG/amg_blas_nse2d.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object src/AMG/CMakeFiles/amg.dir/amg_blas_nse2d.c.o"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/amg.dir/amg_blas_nse2d.c.o   -c /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_blas_nse2d.c

src/AMG/CMakeFiles/amg.dir/amg_blas_nse2d.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/amg.dir/amg_blas_nse2d.c.i"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_blas_nse2d.c > CMakeFiles/amg.dir/amg_blas_nse2d.c.i

src/AMG/CMakeFiles/amg.dir/amg_blas_nse2d.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/amg.dir/amg_blas_nse2d.c.s"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_blas_nse2d.c -o CMakeFiles/amg.dir/amg_blas_nse2d.c.s

src/AMG/CMakeFiles/amg.dir/amg_blas_nse3d.c.o: src/AMG/CMakeFiles/amg.dir/flags.make
src/AMG/CMakeFiles/amg.dir/amg_blas_nse3d.c.o: ../src/AMG/amg_blas_nse3d.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object src/AMG/CMakeFiles/amg.dir/amg_blas_nse3d.c.o"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/amg.dir/amg_blas_nse3d.c.o   -c /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_blas_nse3d.c

src/AMG/CMakeFiles/amg.dir/amg_blas_nse3d.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/amg.dir/amg_blas_nse3d.c.i"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_blas_nse3d.c > CMakeFiles/amg.dir/amg_blas_nse3d.c.i

src/AMG/CMakeFiles/amg.dir/amg_blas_nse3d.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/amg.dir/amg_blas_nse3d.c.s"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_blas_nse3d.c -o CMakeFiles/amg.dir/amg_blas_nse3d.c.s

src/AMG/CMakeFiles/amg.dir/amg_coarsen.c.o: src/AMG/CMakeFiles/amg.dir/flags.make
src/AMG/CMakeFiles/amg.dir/amg_coarsen.c.o: ../src/AMG/amg_coarsen.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object src/AMG/CMakeFiles/amg.dir/amg_coarsen.c.o"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/amg.dir/amg_coarsen.c.o   -c /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_coarsen.c

src/AMG/CMakeFiles/amg.dir/amg_coarsen.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/amg.dir/amg_coarsen.c.i"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_coarsen.c > CMakeFiles/amg.dir/amg_coarsen.c.i

src/AMG/CMakeFiles/amg.dir/amg_coarsen.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/amg.dir/amg_coarsen.c.s"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_coarsen.c -o CMakeFiles/amg.dir/amg_coarsen.c.s

src/AMG/CMakeFiles/amg.dir/amg_iter.c.o: src/AMG/CMakeFiles/amg.dir/flags.make
src/AMG/CMakeFiles/amg.dir/amg_iter.c.o: ../src/AMG/amg_iter.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object src/AMG/CMakeFiles/amg.dir/amg_iter.c.o"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/amg.dir/amg_iter.c.o   -c /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_iter.c

src/AMG/CMakeFiles/amg.dir/amg_iter.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/amg.dir/amg_iter.c.i"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_iter.c > CMakeFiles/amg.dir/amg_iter.c.i

src/AMG/CMakeFiles/amg.dir/amg_iter.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/amg.dir/amg_iter.c.s"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_iter.c -o CMakeFiles/amg.dir/amg_iter.c.s

src/AMG/CMakeFiles/amg.dir/amg_low.c.o: src/AMG/CMakeFiles/amg.dir/flags.make
src/AMG/CMakeFiles/amg.dir/amg_low.c.o: ../src/AMG/amg_low.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object src/AMG/CMakeFiles/amg.dir/amg_low.c.o"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/amg.dir/amg_low.c.o   -c /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_low.c

src/AMG/CMakeFiles/amg.dir/amg_low.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/amg.dir/amg_low.c.i"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_low.c > CMakeFiles/amg.dir/amg_low.c.i

src/AMG/CMakeFiles/amg.dir/amg_low.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/amg.dir/amg_low.c.s"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_low.c -o CMakeFiles/amg.dir/amg_low.c.s

src/AMG/CMakeFiles/amg.dir/amg_solve_main.c.o: src/AMG/CMakeFiles/amg.dir/flags.make
src/AMG/CMakeFiles/amg.dir/amg_solve_main.c.o: ../src/AMG/amg_solve_main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object src/AMG/CMakeFiles/amg.dir/amg_solve_main.c.o"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/amg.dir/amg_solve_main.c.o   -c /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_solve_main.c

src/AMG/CMakeFiles/amg.dir/amg_solve_main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/amg.dir/amg_solve_main.c.i"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_solve_main.c > CMakeFiles/amg.dir/amg_solve_main.c.i

src/AMG/CMakeFiles/amg.dir/amg_solve_main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/amg.dir/amg_solve_main.c.s"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_solve_main.c -o CMakeFiles/amg.dir/amg_solve_main.c.s

src/AMG/CMakeFiles/amg.dir/amg_solvers.c.o: src/AMG/CMakeFiles/amg.dir/flags.make
src/AMG/CMakeFiles/amg.dir/amg_solvers.c.o: ../src/AMG/amg_solvers.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object src/AMG/CMakeFiles/amg.dir/amg_solvers.c.o"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/amg.dir/amg_solvers.c.o   -c /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_solvers.c

src/AMG/CMakeFiles/amg.dir/amg_solvers.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/amg.dir/amg_solvers.c.i"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_solvers.c > CMakeFiles/amg.dir/amg_solvers.c.i

src/AMG/CMakeFiles/amg.dir/amg_solvers.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/amg.dir/amg_solvers.c.s"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_solvers.c -o CMakeFiles/amg.dir/amg_solvers.c.s

src/AMG/CMakeFiles/amg.dir/amg_sp.c.o: src/AMG/CMakeFiles/amg.dir/flags.make
src/AMG/CMakeFiles/amg.dir/amg_sp.c.o: ../src/AMG/amg_sp.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object src/AMG/CMakeFiles/amg.dir/amg_sp.c.o"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/amg.dir/amg_sp.c.o   -c /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_sp.c

src/AMG/CMakeFiles/amg.dir/amg_sp.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/amg.dir/amg_sp.c.i"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_sp.c > CMakeFiles/amg.dir/amg_sp.c.i

src/AMG/CMakeFiles/amg.dir/amg_sp.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/amg.dir/amg_sp.c.s"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_sp.c -o CMakeFiles/amg.dir/amg_sp.c.s

src/AMG/CMakeFiles/amg.dir/amg_variables.c.o: src/AMG/CMakeFiles/amg.dir/flags.make
src/AMG/CMakeFiles/amg.dir/amg_variables.c.o: ../src/AMG/amg_variables.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building C object src/AMG/CMakeFiles/amg.dir/amg_variables.c.o"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/amg.dir/amg_variables.c.o   -c /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_variables.c

src/AMG/CMakeFiles/amg.dir/amg_variables.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/amg.dir/amg_variables.c.i"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_variables.c > CMakeFiles/amg.dir/amg_variables.c.i

src/AMG/CMakeFiles/amg.dir/amg_variables.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/amg.dir/amg_variables.c.s"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG/amg_variables.c -o CMakeFiles/amg.dir/amg_variables.c.s

# Object files for target amg
amg_OBJECTS = \
"CMakeFiles/amg.dir/amg_1d_prec.c.o" \
"CMakeFiles/amg.dir/amg_2d_prec.c.o" \
"CMakeFiles/amg.dir/amg_blas.c.o" \
"CMakeFiles/amg.dir/amg_blas_nse2d.c.o" \
"CMakeFiles/amg.dir/amg_blas_nse3d.c.o" \
"CMakeFiles/amg.dir/amg_coarsen.c.o" \
"CMakeFiles/amg.dir/amg_iter.c.o" \
"CMakeFiles/amg.dir/amg_low.c.o" \
"CMakeFiles/amg.dir/amg_solve_main.c.o" \
"CMakeFiles/amg.dir/amg_solvers.c.o" \
"CMakeFiles/amg.dir/amg_sp.c.o" \
"CMakeFiles/amg.dir/amg_variables.c.o"

# External object files for target amg
amg_EXTERNAL_OBJECTS =

/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/amg_1d_prec.c.o
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/amg_2d_prec.c.o
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/amg_blas.c.o
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/amg_blas_nse2d.c.o
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/amg_blas_nse3d.c.o
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/amg_coarsen.c.o
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/amg_iter.c.o
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/amg_low.c.o
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/amg_solve_main.c.o
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/amg_solvers.c.o
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/amg_sp.c.o
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/amg_variables.c.o
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/build.make
/home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a: src/AMG/CMakeFiles/amg.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Linking C static library /home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a"
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && $(CMAKE_COMMAND) -P CMakeFiles/amg.dir/cmake_clean_target.cmake
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/amg.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/AMG/CMakeFiles/amg.dir/build: /home/divij/Documents/ParMooN_Codes/ParMooN_Output/TBE2D/lib/libamg.a

.PHONY : src/AMG/CMakeFiles/amg.dir/build

src/AMG/CMakeFiles/amg.dir/clean:
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG && $(CMAKE_COMMAND) -P CMakeFiles/amg.dir/cmake_clean.cmake
.PHONY : src/AMG/CMakeFiles/amg.dir/clean

src/AMG/CMakeFiles/amg.dir/depend:
	cd /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/divij/Documents/ParMooN_Codes/ParMooN_CMG /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/src/AMG /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG /home/divij/Documents/ParMooN_Codes/ParMooN_CMG/BUILD/src/AMG/CMakeFiles/amg.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/AMG/CMakeFiles/amg.dir/depend

