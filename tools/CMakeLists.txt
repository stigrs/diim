cmake_minimum_required(VERSION 3.13)
project(DIIM_tools CXX)

include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
link_directories(${BLAS_LIBRARY_DIRS})

if(NOT TARGET scilib::scilib)
    find_package(scilib REQUIRED)
endif()

if(NOT TARGET nlohmann_json::nlohmann_json)
    find_package(nlohmann_json REQUIRED)
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
        nlohmann_json::nlohmann_json
	    ${BLAS_LIBRARIES}
    ) 
    target_compile_options(
        ${name}
        PRIVATE
        "$<$<CONFIG:Debug>:${Diim_CXX_FLAGS_DEBUG}>"
        "$<$<CONFIG:Release>:${Diim_CXX_FLAGS_RELEASE}>"
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
