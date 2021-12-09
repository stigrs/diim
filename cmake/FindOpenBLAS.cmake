# Copyright (c) 2018 Stig Rune Sellevag
#
# This file is distributed under the MIT License. See the accompanying file
# LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
# and conditions.

# Find the OpenBLAS Library.
#
# The following variables are defined if MKL is found:
#
# BLAS_FOUND           - System has OpenBLAS
# BLAS_BINARY_DIRS     - OpenBLAS binary files directory
# BLAS_INCLUDE_DIRS    - OpenBLAS include file directory
# LAPACKE_INCLUDE_DIRS - LAPACKE include file directory
# BLAS_LIBRARY_DIRS    - OpenBLAS library file directory
# BLAS_LIBRARIES       - The OpenBLAS libraries
#
# The environmental variable BLAS_ROOT is used to find the library.
#
# Example usage:
#
# find_package(OpenBLAS)
# if(BLAS_FOUND)
#	  include_directories(${BLAS_INCLUDE_DIRS})
#	  include_directories(${LAPACKE_INCLUDE_DIRS})
#	  link_directories(${BLAS_LIBRARY_DIRS})
#     target_link_libraries(TARGET ${BLAS_LIBRARIES})
# endif()
#
# Note: 
# - Currently, the Intel LP64 interface layer is used for Intel(R) 64
#   architecture.

if(DEFINED ENV{BLAS_ROOT})
    set(BLAS_DIR $ENV{BLAS_ROOT})
else()
    if(WIN32)
        set(BLAS_DIR C:/local/OpenBLAS.0.2.14.1/lib/native)
    elseif(APPLE)
        set(BLAS_DIR /usr/local/opt/openblas)
    else()
        set(BLAS_DIR /usr)
    endif()
endif()
if(WIN32)
	if(${CMAKE_CL_64} EQUAL 1)
        find_file(BLAS_LIBRARIES libopenblas.dll.a HINTS ${BLAS_DIR}/lib/x64)
        set(BLAS_LIBRARY_DIRS "${BLAS_DIR}/lib/x64")
        set(BLAS_BINARY_DIRS "${BLAS_DIR}/bin/x64")
    else()
        find_file(BLAS_LIBRARIES libopenblas.dll.a HINTS ${BLAS_DIR}/lib/win32)
        set(BLAS_LIBRARY_DIRS "${BLAS_DIR}/lib/win32")
        set(BLAS_BINARY_DIRS "${BLAS_DIR}/bin/win32")
    endif()
elseif(APPLE)
    find_file(BLAS_LIBRARIES libopenblas.dylib HINTS ${BLAS_DIR}/lib)
    set(BLAS_LIBRARY_DIRS "${BLAS_DIR}/lib")
    set(BLAS_BINARY_DIRS "${BLAS_DIR}/lib")
else()
    find_file(LIBOPENBLAS libopenblas.so HINTS "${BLAS_DIR}/lib/x86_64-linux-gnu")
    find_file(LIBLAPACKE liblapacke.so HINTS "${BLAS_DIR}/lib/x86_64-linux-gnu")
    set(BLAS_LIBRARIES ${LIBOPENBLAS} ${LIBLAPACKE})
    set(BLAS_LIBRARY_DIRS "${BLAS_DIR}/lib/x86_64-linux-gnu")
    set(BLAS_BINARY_DIRS "${BLAS_DIR}/lib/x86_64-linux-gnu")
endif()

find_path(BLAS_INCLUDE_DIRS cblas.h HINTS ${BLAS_DIR}/include ${BLAS_DIR}/include/x86_64/linux-gnu)
find_path(LAPACKE_INCLUDE_DIRS lapacke.h HINTS ${BLAS_DIR}/include)

if(BLAS_INCLUDE_DIRS AND LAPACKE_INCLUDE_DIRS AND BLAS_LIBRARIES)
    set(BLAS_FOUND ON)
endif()

