# simple CMake file for quicksort experiments

set(SOURCES main.cpp)

set(LIBS rt)

# if have PAPI, set includes, library and -D define switch
if(PAPI_FOUND)
  include_directories(${PAPI_INCLUDE_DIRS})
  add_definitions("-DWITH_PAPI")
  set(LIBS ${LIBS} ${PAPI_LIBRARIES})
endif()

add_executable(hashingtest ${SOURCES})
target_link_libraries(hashingtest ${LIBS})

