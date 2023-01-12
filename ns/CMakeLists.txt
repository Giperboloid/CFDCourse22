cmake_minimum_required(VERSION 3.0)

Project(cfdapp)

set(CMAKE_CXX_STANDARD 17)

# build all binaries to build/bin
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR}/lib)
include_directories(${BOOST_ROOT})
include_directories(${AMGCL_ROOT})

set(CFDLIB_TARGET cfdlib)
#
#find_package(OpenMP)
#if (OPENMP_FOUND)
#	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

# default input directory to "input_data/"
add_definitions(-DCFDLIB_INPUT_DIR_PATH="${CMAKE_SOURCE_DIR}/input_data/")

# default output directory
if (MSVC)
	# build/output for visual studio
	set(OUTPUT_DIR_PATH ${CMAKE_BINARY_DIR}/output/)
	add_definitions(-DCFDLIB_OUTPUT_DIR_PATH="${OUTPUT_DIR_PATH}")
	file(MAKE_DIRECTORY ${OUTPUT_DIR_PATH})
else(MSVC)
	# ./ for linux
	add_definitions(-DCFDLIB_OUTPUT_DIR_PATH="")
endif(MSVC)

add_subdirectory(lib)
add_subdirectory(progs)
