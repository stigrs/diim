cmake_minimum_required(VERSION 3.4.3)
project(libdiim CXX)

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

set(
    SRC_FILES
    diim.cpp
)
add_library(
    diim 
    ${SRC_FILES}
)
add_dependencies(diim Stdutils Numlib)
target_link_libraries(
    diim 
    ${BLAS_LIBRARIES}
    ${Numlib_LIBRARIES} 
) 

install(
    TARGETS diim
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

