#################################
# cmake file for building DATABASE
# @author LAGARDE François
#################################

# project name
PROJECT( Streamout )

# project version
SET( STREAMOUT_VERSION_MAJOR 1 )
SET( STREAMOUT_VERSION_MINOR 0 )
SET( STREAMOUT_VERSION_PATCH 0 )
file(COPY ./include/ DESTINATION ${PROJECT_INCLUDE_DIR})
ADD_DEFINITIONS(-std=c++1y)
FILE (GLOB_RECURSE STREAMOUT_SRC_FILES ./src/*)
ADD_LIBRARY( Streamout SHARED ${STREAMOUT_SRC_FILES} )
TARGET_LINK_LIBRARIES(Streamout Reader ${ROOT_LIBRARIES})
INSTALL(TARGETS Streamout DESTINATION ${CMAKE_SOURCE_DIR}/lib ) 
