cmake_minimum_required(VERSION 3.17)
project(DIIM_tools CXX)

include_directories(${DIIM_SOURCE_DIR}/include)
include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})

if(NOT TARGET std::mdspan)
  find_package(mdspan REQUIRED)
endif()

if(NOT TARGET Scilib::Scilib)
    find_package(scilib REQUIRED)
endif()

if(NOT TARGET Stdutils::Stdutils)
    find_package(stdutils REQUIRED)
endif()

if(NOT TARGET Diim::Diim)
    find_package(diim REQUIRED)
endif()

function(add_tool name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries (
        ${name} 
        PRIVATE
        Diim::Diim
        Scilib::Scilib
        Stdutils::Stdutils
        std::mdspan
	    ${BLAS_LIBRARIES}
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
