# CMake file for hashing project

project(hashing)

cmake_minimum_required(VERSION 2.8)

option(WITH_PAPI "Use PAPI library for performance counting" ON)

# disallow in-source builds

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")


# Enable warnings

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")

#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -g -funroll-loops")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -g")
#string(REPLACE "-O3" "-O2" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

find_package(Boost COMPONENTS random REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Find PAPI library

if(WITH_PAPI)

  find_path(PAPI_INCLUDE_DIRS papi.h
    PATHS /usr/include/ /usr/local/include/)

  find_library(PAPI_LIBRARIES NAMES papi
    PATHS /usr/lib/ /usr/local/lib/)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(PAPI DEFAULT_MSG
    PAPI_INCLUDE_DIRS
    PAPI_LIBRARIES
    )

endif()

# descend into source
add_subdirectory(src)
