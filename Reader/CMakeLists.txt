#################################
# cmake file for building DATABASE
# @author LAGARDE François
#################################

# project name
PROJECT( Reader )

# project version
SET( COMMON_VERSION_MAJOR 1 )
SET( COMMON_VERSION_MINOR 0 )
SET( COMMON_VERSION_PATCH 0 )
ADD_DEFINITIONS(-std=c++11)
file(COPY ./include/ DESTINATION ${PROJECT_INCLUDE_DIR})
FILE (GLOB_RECURSE READER_SRC_FILES ./src/*)
ADD_LIBRARY( Reader STATIC ${READER_SRC_FILES} )
TARGET_LINK_LIBRARIES(Reader Geometry Tinyxml )
INSTALL(TARGETS Reader DESTINATION ${LIBRARY_OUTPUT_PATH} ) 

