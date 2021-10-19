# Copyright (c) 2018 Stig Rune Sellevag
#
# This file is distributed under the MIT License. See the accompanying file
# LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
# and conditions.

# Find the Numlib Library.
#
# The following variables are defined if Numlib is found:
#
# Numlib_FOUND        - System has Numlib
# Numlib_INCLUDE_DIRS - Numlib include file directory
# Numlib_LIBRARY_DIRS - Numlib library file directory
# Numlib_LIBRARIES    - The Numlib libraries
#
# The environmental variable Numlib_ROOT is used to find the library.
#
# Example usage:
#
# find_package(Numlib)
# if(Numlib_FOUND)
#     include_directories(${Numlib_INCLUDE_DIRS})
#     link_directories(${Numlib_LIBRARY_DIRS})
#     target_link_libraries(TARGET ${Numlib_LIBRARIES})
# endif()

find_path(Numlib_INCLUDE_DIRS numlib HINTS $ENV{Numlib_ROOT}/include $ENV{HOME}/include)
if(WIN32)
    find_path(Numlib_LIBRARY_DIRS num.lib HINTS $ENV{Numlib_ROOT}/lib $ENV{HOME}/lib)
else()
    if(APPLE)
        find_path(Numlib_LIBRARY_DIRS libnum.dylib HINTS $ENV{Numlib_ROOT}/lib $ENV{HOME}/lib)
    else()
	find_path(Numlib_LIBRARY_DIRS libnum.so HINTS $ENV{Numlib_ROOT}/lib $ENV{HOME}/lib)
    endif()
endif()

if(Numlib_INCLUDE_DIRS AND Numlib_LIBRARY_DIRS)
    set(Numlib_LIBRARIES "num")
    set(Numlib_FOUND ON)
    message("-- Numlib found")
endif()
