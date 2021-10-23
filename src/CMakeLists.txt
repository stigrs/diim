cmake_minimum_required(VERSION 3.4.3)
project(DIIM_apps CXX)

include_directories(${DIIM_SOURCE_DIR}/include)
include_directories(${Stdutils_INCLUDE_DIRS})
include_directories(${Numlib_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
link_directories(${DIIM_BINARY_DIR})
link_directories(${Numlib_LIBRARY_DIRS})
if(MKL_FOUND)
    link_directories(${MKL_LIBRARY_DIRS})
endif()

if(GIT_FOUND)
    # Add cxxopts
    list(APPEND CXXOPTS_CMAKE_ARGS
        "-DCMAKE_INSTALL_PREFIX:STRING=${CMAKE_BINARY_DIR}/external"
        "-DCXXOPTS_BUILD_EXAMPLES=OFF"
        "-DCXXOPTS_BUILD_TESTS=OFF"
    )
    ExternalProject_Add(
        cxxopts
        PREFIX ${CMAKE_BINARY_DIR}/cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts
        GIT_TAG v2.1.1
        CMAKE_ARGS ${CXXOPTS_CMAKE_ARGS}
        LOG_DOWNLOAD 1
        UPDATE_DISCONNECTED 1
    )
else()
    # Assume catch is installed in a system directory
    add_custom_target(cxxopts)
endif()
 
# Prepare cxxopts library for other executables.
add_library(Cxxopts INTERFACE)
target_include_directories(
    Cxxopts INTERFACE 
    ${CMAKE_BINARY_DIR}/external/include
)

function(add_app name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries (
        ${name} 
        iim
	    ${BLAS_LIBRARIES}
        ${Numlib_LIBRARIES}
        Cxxopts
    ) 
    add_dependencies(${name} cxxopts Stdutils Numlib)
endfunction()

set(PROGRAMS 
    diim
)

foreach(program ${PROGRAMS})
    add_app(${program})
endforeach()

install(
    TARGETS ${PROGRAMS}
    DESTINATION bin  
)

