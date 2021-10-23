cmake_minimum_required(VERSION 3.4.3)
project(DIIM_tests CXX)

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

# Copy *.inp files to build directory.
file(GLOB INP_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.inp")
foreach(inp_file_name ${INP_FILES})
    file(COPY "${inp_file_name}" DESTINATION "${CMAKE_BINARY_DIR}/tests")
endforeach()

# Copy *.csv files to build directory.
file(GLOB CSV_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.csv")
foreach(csv_file_name ${CSV_FILES})
    file(COPY "${csv_file_name}" DESTINATION "${CMAKE_BINARY_DIR}/tests")
endforeach()

if(GIT_FOUND)
    # Add catch
    list(APPEND CATCH_CMAKE_ARGS
        "-DCMAKE_INSTALL_PREFIX:STRING=${CMAKE_BINARY_DIR}/external"
        "-DCATCH_BUILD_TESTING=OFF"
    )
    ExternalProject_Add(
        catch2
        PREFIX ${CMAKE_BINARY_DIR}/catch
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v2.3.0
        CMAKE_ARGS ${CATCH_CMAKE_ARGS}
        LOG_DOWNLOAD 1
        UPDATE_DISCONNECTED 1
    )
else()
    # Assume catch is installed in a system directory
    add_custom_target(catch2)
endif()

# Prepare Catch library for other executables.
add_library(Catch2 INTERFACE)
target_include_directories(
    Catch2 INTERFACE 
    ${CMAKE_BINARY_DIR}/external/include
)

# Create the main executable for each test.
add_library(test_catch2 STATIC test.cpp)
target_link_libraries(
    test_catch2
    Catch2
)
add_dependencies(test_catch2 catch2)

function(add_diim_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries (
        ${name} 
        test_catch2
        iim
	${BLAS_LIBRARIES}
        ${Numlib_LIBRARIES}
    ) 
    add_dependencies(${name} Catch2 Stdutils Numlib)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

set(PROGRAMS 
    test_diim
)

foreach(program ${PROGRAMS})
    add_diim_test(${program})
	if(WIN32)
	    if(NOT MKL_FOUND)
            set_tests_properties(${program} PROPERTIES ENVIRONMENT "PATH=${BLAS_BINARY_DIRS};$ENV{PATH}")
	    endif()
    else()
	    if(NOT MKL_FOUND)
            set_tests_properties(${program} PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${BLAS_BINARY_DIRS}")
	    endif()
	endif()
endforeach()
