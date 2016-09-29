#################################
# cmake file for building DATABASE
# @author LAGARDE François
#################################

# project name
PROJECT( tinyxml )

# project version
SET( TINYXML_VERSION_MAJOR 1 )
SET( TINYXML_VERSION_MINOR 0 )
SET( TINYXML_VERSION_PATCH 0 )
file(COPY ./include/ DESTINATION ${PROJECT_INCLUDE_DIR})
ADD_DEFINITIONS(-std=c++1y)
FILE (GLOB_RECURSE TINYXML_SRC_FILES ./src/*)
ADD_LIBRARY( Tinyxml STATIC ${TINYXML_SRC_FILES} )
TARGET_LINK_LIBRARIES(Tinyxml Reader ${ROOT_LIBRARIES})
INSTALL(TARGETS Tinyxml DESTINATION ${CMAKE_SOURCE_DIR}/lib ) 
