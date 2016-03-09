########################################################
# cmake file for building EventBuilder
# @author Lagarde François, IPNL
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
########################################################

# project name
project(EventBuilder)
include(CheckIncludeFiles)
include(CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  else()
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX11)
	  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	  else()
    if(COMPILER_SUPPORTS_CXX0X)
      CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
      else()
      message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
    endif()
  endif()
endif()


add_definitions(-DTINYXML_USE_STL)

find_library(LOCCI_SUPPORT NAMES occi)
if(LOCCI_SUPPORT)
  message("LOCCI found")
  SET(HaveOCCI 1)
  SET(GCC_COVERAGE_LINK_FLAGS "-locci -lclntsh -g0")
  add_definitions(-DOCCI_SUPPORT)
else()
  message("No LOCCI support")
  SET(HaveOCCI 0)
endif()



SET (CMAKE_SHARED_LINKER_FLAGS "${GCC_COVERAGE_LINK_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${GCC_COVERAGE_LINK_FLAGS}" )
include_directories(include)

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 3 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH a )

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

FIND_PACKAGE( ROOT REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )
ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )


#find_package(ORACLE REQUIRED)
#include_directories(${ORACLE_INCLUDE_DIRS})
#LINK_LIBRARIES( ${ORACLE_LIBRARIES} )
#ADD_DEFINITIONS( ${ORACLE_DEFINITIONS} )


### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF)

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()

### SUB DIRECTORIES #########################################################
ADD_SUBDIRECTORY( source )
add_executable(Plot ./source/src/Plot/Plot.cpp )
# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

