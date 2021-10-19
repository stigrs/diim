# Copyright (c) 2018 Stig Rune Sellevag
#
# This file is distributed under the MIT License. See the accompanying file
# LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
# and conditions.

# Find the Stdutils Library.
#
# The following variables are defined if Stdutils is found:
#
# Stdutils_FOUND        - System has Stdutils
# Stdutils_INCLUDE_DIRS - Stdutils include file directory
#
# The environmental variable Stdutils_ROOT is used to find the library.
#
# Example usage:
#
# find_package(Stdutils)
# if(Stdutils_FOUND)
#	  include_directories(${Stdutils_INCLUDE_DIRS})
# endif()
#

find_path(Stdutils_INCLUDE_DIRS stdutils HINTS $ENV{Stdutils_ROOT}/include $ENV{HOME}/include)
if(Stdutils_INCLUDE_DIRS)
    set(Stdutils_FOUND ON)
	message("-- Stdutils found")
endif()
