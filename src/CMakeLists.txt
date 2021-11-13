cmake_minimum_required(VERSION 3.4.3)
project(DIIM_apps CXX)

include_directories(${DIIM_SOURCE_DIR}/include)
include_directories(${Stdutils_INCLUDE_DIRS})
include_directories(${Numlib_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
link_directories(${DIIM_BINARY_DIR})
link_directories(${Numlib_LIBRARY_DIRS})
if(MKL_FOUND)
    link_directories(${MKL_LIBRARY_DIRS})
endif()

function(add_app name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries (
        ${name} 
        iim
	    ${BLAS_LIBRARIES}
        ${Numlib_LIBRARIES}
    ) 
    add_dependencies(${name} Stdutils Numlib)
endfunction()

set(PROGRAMS 
    diim
    diim_edit
)

foreach(program ${PROGRAMS})
    add_app(${program})
endforeach()

install(
    TARGETS ${PROGRAMS}
    DESTINATION bin  
)
