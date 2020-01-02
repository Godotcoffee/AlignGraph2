cmake_minimum_required(VERSION 3.1.0)
project(PAGraph)

set(CMAKE_CXX_STANDARD 11)

set(LIB_COMPILE_WARNINGS "-Wall -Wno-unused-function -Wno-format")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(debug 1)
endif()

if(debug)
    #set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS} -g -p -pg")
    set(LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS} -g -p ${LIB_COMPILE_WARNINGS}")
    set(CMAKE_BUILD_TYPE Debug) # else CMake adds DNDEBUG
    message("-- COMPILATION IN DEBUG MODE")
else()
    set(LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS} -O3 -DNDEBUG ${LIB_COMPILE_WARNINGS}")
endif()

set(project-libraries pthread)

set(THIRDPARTY_SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty)
include_directories(${THIRDPARTY_SOURCE_DIR})

set(TOOLS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/tools)
include_directories(${TOOLS_SOURCE_DIR})
file(GLOB_RECURSE ToolsFiles ${TOOLS_SOURCE_DIR}/*.cpp)

set(MAIN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/main)
file(GLOB TargetFiles ${MAIN_SOURCE_DIR}/*.cpp)

FOREACH(target_file ${TargetFiles})
    get_filename_component(target ${target_file} NAME_WE)
    add_executable(${target} ${target_file} ${ToolsFiles})
    target_link_libraries(${target} ${project-libraries})
ENDFOREACH()
