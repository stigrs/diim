cmake_minimum_required(VERSION 3.14)
project(libdiim CXX)

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

set(SRC_FILES
    diim.cpp
    perturbation.cpp
    utils.cpp
)
add_library(
    diim 
    ${SRC_FILES}
)
target_link_libraries(
    diim 
    ${BLAS_LIBRARIES}
    ${Scilib_LIBRARIES} 
    std::mdspan
) 

install(
    TARGETS diim
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
