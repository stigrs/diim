cmake_minimum_required(VERSION 3.14)
project(DIIM_tools CXX)

include_directories(${DIIM_SOURCE_DIR}/include)
include_directories(${Stdutils_INCLUDE_DIRS})
include_directories(${Scilib_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
link_directories(${DIIM_BINARY_DIR})
link_directories(${BLAS_LIBRARY_DIRS})

if(NOT TARGET std::mdspan)
  find_package(mdspan REQUIRED)
endif()

function(add_tool name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries (
        ${name} 
        diim
	    ${BLAS_LIBRARIES}
        ${Scilib_LIBRARIES}
        std::mdspan
    ) 
endfunction()

set(PROGRAMS 
    diim_run
    diim_edit
)

foreach(program ${PROGRAMS})
    add_tool(${program})
endforeach()

install(
    TARGETS ${PROGRAMS}
    DESTINATION bin  
)
