# CMAKE generated file: DO NOT EDIT!
# Generated by "NMake Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

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

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE
NULL=nul
!ENDIF
SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2023.2.2\bin\cmake\win\x64\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2023.2.2\bin\cmake\win\x64\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

# Include any dependencies generated for this target.
include easy3d\util\CMakeFiles\easy3d_util.dir\depend.make
# Include any dependencies generated by the compiler for this target.
include easy3d\util\CMakeFiles\easy3d_util.dir\compiler_depend.make

# Include the progress variables for this target.
include easy3d\util\CMakeFiles\easy3d_util.dir\progress.make

# Include the compile flags for this target's objects.
include easy3d\util\CMakeFiles\easy3d_util.dir\flags.make

easy3d\util\CMakeFiles\easy3d_util.dir\chrono_watch.cpp.obj: easy3d\util\CMakeFiles\easy3d_util.dir\flags.make
easy3d\util\CMakeFiles\easy3d_util.dir\chrono_watch.cpp.obj: C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\chrono_watch.cpp
easy3d\util\CMakeFiles\easy3d_util.dir\chrono_watch.cpp.obj: easy3d\util\CMakeFiles\easy3d_util.dir\compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object easy3d/util/CMakeFiles/easy3d_util.dir/chrono_watch.cpp.obj"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	$(CMAKE_COMMAND) -E cmake_cl_compile_depends --dep-file=CMakeFiles\easy3d_util.dir\chrono_watch.cpp.obj.d --working-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util --filter-prefix="注意: 包含文件:  " -- C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /showIncludes /FoCMakeFiles\easy3d_util.dir\chrono_watch.cpp.obj /FdCMakeFiles\easy3d_util.dir\easy3d_util.pdb /FS -c C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\chrono_watch.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\chrono_watch.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_util.dir/chrono_watch.cpp.i"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe > CMakeFiles\easy3d_util.dir\chrono_watch.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\chrono_watch.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\chrono_watch.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_util.dir/chrono_watch.cpp.s"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_util.dir\chrono_watch.cpp.s /c C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\chrono_watch.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\dialogs.cpp.obj: easy3d\util\CMakeFiles\easy3d_util.dir\flags.make
easy3d\util\CMakeFiles\easy3d_util.dir\dialogs.cpp.obj: C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\dialogs.cpp
easy3d\util\CMakeFiles\easy3d_util.dir\dialogs.cpp.obj: easy3d\util\CMakeFiles\easy3d_util.dir\compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object easy3d/util/CMakeFiles/easy3d_util.dir/dialogs.cpp.obj"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	$(CMAKE_COMMAND) -E cmake_cl_compile_depends --dep-file=CMakeFiles\easy3d_util.dir\dialogs.cpp.obj.d --working-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util --filter-prefix="注意: 包含文件:  " -- C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /showIncludes /FoCMakeFiles\easy3d_util.dir\dialogs.cpp.obj /FdCMakeFiles\easy3d_util.dir\easy3d_util.pdb /FS -c C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\dialogs.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\dialogs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_util.dir/dialogs.cpp.i"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe > CMakeFiles\easy3d_util.dir\dialogs.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\dialogs.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\dialogs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_util.dir/dialogs.cpp.s"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_util.dir\dialogs.cpp.s /c C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\dialogs.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\file_system.cpp.obj: easy3d\util\CMakeFiles\easy3d_util.dir\flags.make
easy3d\util\CMakeFiles\easy3d_util.dir\file_system.cpp.obj: C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\file_system.cpp
easy3d\util\CMakeFiles\easy3d_util.dir\file_system.cpp.obj: easy3d\util\CMakeFiles\easy3d_util.dir\compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object easy3d/util/CMakeFiles/easy3d_util.dir/file_system.cpp.obj"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	$(CMAKE_COMMAND) -E cmake_cl_compile_depends --dep-file=CMakeFiles\easy3d_util.dir\file_system.cpp.obj.d --working-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util --filter-prefix="注意: 包含文件:  " -- C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /showIncludes /FoCMakeFiles\easy3d_util.dir\file_system.cpp.obj /FdCMakeFiles\easy3d_util.dir\easy3d_util.pdb /FS -c C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\file_system.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\file_system.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_util.dir/file_system.cpp.i"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe > CMakeFiles\easy3d_util.dir\file_system.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\file_system.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\file_system.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_util.dir/file_system.cpp.s"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_util.dir\file_system.cpp.s /c C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\file_system.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\stop_watch.cpp.obj: easy3d\util\CMakeFiles\easy3d_util.dir\flags.make
easy3d\util\CMakeFiles\easy3d_util.dir\stop_watch.cpp.obj: C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\stop_watch.cpp
easy3d\util\CMakeFiles\easy3d_util.dir\stop_watch.cpp.obj: easy3d\util\CMakeFiles\easy3d_util.dir\compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object easy3d/util/CMakeFiles/easy3d_util.dir/stop_watch.cpp.obj"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	$(CMAKE_COMMAND) -E cmake_cl_compile_depends --dep-file=CMakeFiles\easy3d_util.dir\stop_watch.cpp.obj.d --working-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util --filter-prefix="注意: 包含文件:  " -- C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /showIncludes /FoCMakeFiles\easy3d_util.dir\stop_watch.cpp.obj /FdCMakeFiles\easy3d_util.dir\easy3d_util.pdb /FS -c C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\stop_watch.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\stop_watch.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_util.dir/stop_watch.cpp.i"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe > CMakeFiles\easy3d_util.dir\stop_watch.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\stop_watch.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\stop_watch.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_util.dir/stop_watch.cpp.s"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_util.dir\stop_watch.cpp.s /c C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\stop_watch.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\string.cpp.obj: easy3d\util\CMakeFiles\easy3d_util.dir\flags.make
easy3d\util\CMakeFiles\easy3d_util.dir\string.cpp.obj: C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\string.cpp
easy3d\util\CMakeFiles\easy3d_util.dir\string.cpp.obj: easy3d\util\CMakeFiles\easy3d_util.dir\compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object easy3d/util/CMakeFiles/easy3d_util.dir/string.cpp.obj"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	$(CMAKE_COMMAND) -E cmake_cl_compile_depends --dep-file=CMakeFiles\easy3d_util.dir\string.cpp.obj.d --working-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util --filter-prefix="注意: 包含文件:  " -- C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /showIncludes /FoCMakeFiles\easy3d_util.dir\string.cpp.obj /FdCMakeFiles\easy3d_util.dir\easy3d_util.pdb /FS -c C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\string.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\string.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_util.dir/string.cpp.i"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe > CMakeFiles\easy3d_util.dir\string.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\string.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\string.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_util.dir/string.cpp.s"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_util.dir\string.cpp.s /c C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\string.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\threading.cpp.obj: easy3d\util\CMakeFiles\easy3d_util.dir\flags.make
easy3d\util\CMakeFiles\easy3d_util.dir\threading.cpp.obj: C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\threading.cpp
easy3d\util\CMakeFiles\easy3d_util.dir\threading.cpp.obj: easy3d\util\CMakeFiles\easy3d_util.dir\compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object easy3d/util/CMakeFiles/easy3d_util.dir/threading.cpp.obj"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	$(CMAKE_COMMAND) -E cmake_cl_compile_depends --dep-file=CMakeFiles\easy3d_util.dir\threading.cpp.obj.d --working-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util --filter-prefix="注意: 包含文件:  " -- C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /showIncludes /FoCMakeFiles\easy3d_util.dir\threading.cpp.obj /FdCMakeFiles\easy3d_util.dir\easy3d_util.pdb /FS -c C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\threading.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\threading.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_util.dir/threading.cpp.i"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe > CMakeFiles\easy3d_util.dir\threading.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\threading.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

easy3d\util\CMakeFiles\easy3d_util.dir\threading.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_util.dir/threading.cpp.s"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_util.dir\threading.cpp.s /c C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util\threading.cpp
<<
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

# Object files for target easy3d_util
easy3d_util_OBJECTS = \
"CMakeFiles\easy3d_util.dir\chrono_watch.cpp.obj" \
"CMakeFiles\easy3d_util.dir\dialogs.cpp.obj" \
"CMakeFiles\easy3d_util.dir\file_system.cpp.obj" \
"CMakeFiles\easy3d_util.dir\stop_watch.cpp.obj" \
"CMakeFiles\easy3d_util.dir\string.cpp.obj" \
"CMakeFiles\easy3d_util.dir\threading.cpp.obj"

# External object files for target easy3d_util
easy3d_util_EXTERNAL_OBJECTS =

lib\easy3d_util.lib: easy3d\util\CMakeFiles\easy3d_util.dir\chrono_watch.cpp.obj
lib\easy3d_util.lib: easy3d\util\CMakeFiles\easy3d_util.dir\dialogs.cpp.obj
lib\easy3d_util.lib: easy3d\util\CMakeFiles\easy3d_util.dir\file_system.cpp.obj
lib\easy3d_util.lib: easy3d\util\CMakeFiles\easy3d_util.dir\stop_watch.cpp.obj
lib\easy3d_util.lib: easy3d\util\CMakeFiles\easy3d_util.dir\string.cpp.obj
lib\easy3d_util.lib: easy3d\util\CMakeFiles\easy3d_util.dir\threading.cpp.obj
lib\easy3d_util.lib: easy3d\util\CMakeFiles\easy3d_util.dir\build.make
lib\easy3d_util.lib: easy3d\util\CMakeFiles\easy3d_util.dir\objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX static library ..\..\lib\easy3d_util.lib"
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	$(CMAKE_COMMAND) -P CMakeFiles\easy3d_util.dir\cmake_clean_target.cmake
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	C:\PROGRA~1\MIB055~1\2022\COMMUN~1\VC\Tools\MSVC\1437~1.328\bin\Hostx64\x64\lib.exe /nologo /machine:x64 /out:..\..\lib\easy3d_util.lib @CMakeFiles\easy3d_util.dir\objects1.rsp 
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug

# Rule to build all files generated by this target.
easy3d\util\CMakeFiles\easy3d_util.dir\build: lib\easy3d_util.lib
.PHONY : easy3d\util\CMakeFiles\easy3d_util.dir\build

easy3d\util\CMakeFiles\easy3d_util.dir\clean:
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util
	$(CMAKE_COMMAND) -P CMakeFiles\easy3d_util.dir\cmake_clean.cmake
	cd C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug
.PHONY : easy3d\util\CMakeFiles\easy3d_util.dir\clean

easy3d\util\CMakeFiles\easy3d_util.dir\depend:
	$(CMAKE_COMMAND) -E cmake_depends "NMake Makefiles" C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\easy3d\util C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util C:\Users\Acer\Documents\GitHub\Geo1016\A1_Calibration_Code\cmake-build-debug\easy3d\util\CMakeFiles\easy3d_util.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : easy3d\util\CMakeFiles\easy3d_util.dir\depend

