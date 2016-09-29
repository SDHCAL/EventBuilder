#################################
# cmake file for building DATABASE
# @author LAGARDE François
#################################

# project name
PROJECT( Geometry )

# project version
SET( GEOMETRY_VERSION_MAJOR 1 )
SET( GEOMETRY_VERSION_MINOR 0 )
SET( GEOMETRY_VERSION_PATCH 0 )

file(COPY ./include/ DESTINATION ${PROJECT_INCLUDE_DIR})
ADD_DEFINITIONS(-std=c++11)
FILE (GLOB_RECURSE GEOMETRY_SRC_FILES ./src/*)
ADD_LIBRARY( Geometry STATIC ${GEOMETRY_SRC_FILES} )
INSTALL(TARGETS Geometry DESTINATION ${LIBRARY_OUTPUT_PATH} ) 
