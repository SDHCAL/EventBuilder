########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT(MyProgramm)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
include_directories(include)

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 9 )
SET( ${PROJECT_NAME}_VERSION_PATCH 9 )



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


# optional package
#FIND_PACKAGE( AIDA )
#IF( AIDA_FOUND )
#    INCLUDE_DIRECTORIES( ${AIDA_INCLUDE_DIRS} )
#    LINK_LIBRARIES( ${AIDA_LIBRARIES} )
#    ADD_DEFINITIONS( "-DMARLIN_USE_AIDA" )
#    MESSAGE( STATUS "AIDA -- found" )
#ELSE()
#    MESSAGE( STATUS "AIDA -- not found" )
#ENDIF()



### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" ON )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



ADD_SUBDIRECTORY( ./Common)
ADD_SUBDIRECTORY( ./StreamOut)
ADD_SUBDIRECTORY( ./Trivent )
ADD_SUBDIRECTORY(./Analysis)
ADD_SUBDIRECTORY(./Noise)


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

