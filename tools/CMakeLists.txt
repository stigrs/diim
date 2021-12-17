cmake_minimum_required(VERSION 3.17)
project(DIIM_tools CXX)

include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
link_directories(${BLAS_LIBRARY_DIRS})

if(NOT TARGET std::mdspan)
  find_package(mdspan REQUIRED)
endif()

if(NOT TARGET scilib::scilib)
    find_package(scilib REQUIRED)
endif()

if(NOT TARGET stdutils::stdutils)
    find_package(stdutils REQUIRED)
endif()

if(NOT TARGET diim::diim)
    find_package(diim REQUIRED)
endif()

function(add_tool name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries (
        ${name} 
        PRIVATE
        diim::diim
        scilib::scilib
        stdutils::stdutils
        std::mdspan
	    ${BLAS_LIBRARIES}
    ) 
endfunction()

set(PROGRAMS 
    diim_edit
    diim_gen
    diim_run
)

foreach(program ${PROGRAMS})
    add_tool(${program})
endforeach()

install(
    TARGETS ${PROGRAMS}
    DESTINATION bin  
)
