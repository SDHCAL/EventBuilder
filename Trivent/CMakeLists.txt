#################################
# cmake file for building DATABASE
# @author LAGARDE François
#################################

# project name
PROJECT( Trivent )

# project version
SET( TRIVENT_VERSION_MAJOR 1 )
SET( TRIVENT_VERSION_MINOR 0 )
SET( TRIVENT_VERSION_PATCH 0 )
file(COPY ./include/ DESTINATION ${PROJECT_INCLUDE_DIR})
ADD_DEFINITIONS(-std=c++1y)
FILE (GLOB_RECURSE TRIVENT_SRC_FILES ./src/*)
ADD_LIBRARY( Trivent SHARED ${TRIVENT_SRC_FILES} )
TARGET_LINK_LIBRARIES(Trivent Reader ${ROOT_LIBRARIES})
INSTALL(TARGETS Trivent DESTINATION ${CMAKE_SOURCE_DIR}/lib ) 
