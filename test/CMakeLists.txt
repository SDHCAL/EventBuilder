########################################################
# cmake file for building Marlin
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
########################################################

# project name
PROJECT( Test )
# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 8 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )
ADD_DEFINITIONS(-std=c++1y)
# include directories
set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME ON) 
find_package(Boost 1.45.0 COMPONENTS filesystem) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    add_executable(Marlin ./src/main.cpp) 
    target_link_libraries(Marlin ${Boost_LIBRARIES} Tinyxml Marlino ${CMAKE_DL_LIBS})
endif()
INSTALL( TARGETS Marlin DESTINATION ${CMAKE_BINARY_DIR} )
# ----------------------------------------------------------------------------
