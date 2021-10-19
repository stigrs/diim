cmake_minimum_required(VERSION 3.4.3)
project(DIIM CXX)

#-------------------------------------------------------------------------------
#
# Set compiler options:

if(UNIX AND NOT APPLE)
    find_package(OpenMP)
endif()

# Enforce C++14 standard.
set(CMAKE_CXX_STANDARD 14)

# Export compile commands.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(WIN32)
    option(BUILD_SHARED_LIBS "Compile shared libraries." OFF)
else()
    option(BUILD_SHARED_LIBS "Compile shared libraries." ON)
endif()

# Set default MSVC compiler options to avoid D9025 error.
if(MSVC)
    set(CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS /D_CRT_SECURE_NO_WARNINGS")
endif()

# Set release compiler options.
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wpedantic -Wall -Wshadow -Wextra -Wno-maybe-uninitialized -Wno-missing-braces -Werror -Wno-error=strict-overflow")
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wpedantic -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/DNDEBUG /std:c++14 /W4 /WX /wd4100 /GR /EHsc /MT /O2")
endif()

# Set debugging compiler options.
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wpedantic -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow")
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wpedantic -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "/std:c++14 /EHsc /W4 /WX /wd4100 /MT /GR /Zi /Od")
endif()

# Set OpenMP compiler options.
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#-------------------------------------------------------------------------------
#
# Include external dependencies:

include(ExternalProject)
find_package(Git)

# Boost is required.
find_package(Boost 1.54 REQUIRED)

# BLAS and LAPACK is required.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
if(DEFINED ENV{MKLROOT})
	find_package(MKL)
	if(MKL_FOUND)
		add_definitions(${MKL_DEFINITIONS})
		set(BLAS_INCLUDE_DIRS ${MKL_INCLUDE_DIRS})
		set(BLAS_LIBRARIES ${MKL_LIBRARIES})
	endif()
else()
    find_package(OpenBLAS)
endif()

# Stdutils is required.
find_package(Stdutils)
if(Stdutils_FOUND)
    add_custom_target(Stdutils)
else()
    list(APPEND STDUTILS_CMAKE_ARGS
        "-DCMAKE_INSTALL_PREFIX:STRING=${CMAKE_BINARY_DIR}/external"
        "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
        "-DBUILD_TESTS=OFF"
    )
    ExternalProject_Add(
        Stdutils
        PREFIX ${CMAKE_BINARY_DIR}/stdutils
        GIT_REPOSITORY https://github.com/stigrs/stdutils.git
        CMAKE_ARGS ${STDUTILS_CMAKE_ARGS}
        LOG_DOWNLOAD 1
        UPDATE_DISCONNECTED 0
    )
    set(Stdutils_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/external/include)
endif()

# Numlib is required.
find_package(Numlib)
if(Numlib_FOUND)
    add_custom_target(Numlib)
else()
    list(APPEND NUMLIB_CMAKE_ARGS
        "-DCMAKE_INSTALL_PREFIX:STRING=${CMAKE_BINARY_DIR}/external"
        "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
        "-DBUILD_TESTS=OFF"
        "-DBUILD_BENCH=OFF"
        "-DBUILD_EXAMPLES=OFF"
    )
    ExternalProject_Add(
        Numlib
        PREFIX ${CMAKE_BINARY_DIR}/numlib
        GIT_REPOSITORY https://github.com/stigrs/numlib.git
        CMAKE_ARGS ${NUMLIB_CMAKE_ARGS}
        LOG_DOWNLOAD 1
        UPDATE_DISCONNECTED 0
    )
    set(Numlib_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/external/include)
    set(Numlib_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/external/lib)
    set(Numlib_LIBRARIES "num")
endif()

#-------------------------------------------------------------------------------
#
# Build the software:

add_subdirectory(lib)
add_subdirectory(src)

# Enable testing.
option(BUILD_TESTS "Build tests." ON)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Install target.
install(
    DIRECTORY include/diim
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)
