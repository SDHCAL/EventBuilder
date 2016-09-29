#################################
# cmake file for building DATABASE
# @author LAGARDE François
#################################

# project name
PROJECT( Database )

# project version
SET( DATABASE_VERSION_MAJOR 1 )
SET( DATABASE_VERSION_MINOR 0 )
SET( DATABASE_VERSION_PATCH 0 )

file(COPY ./include/ DESTINATION ${PROJECT_INCLUDE_DIR})
ADD_DEFINITIONS(-std=c++11)
FILE (GLOB_RECURSE DATABASE_SRC_FILES ./src/*)
ADD_LIBRARY( Database STATIC ${DATABASE_SRC_FILES} )
INSTALL(TARGETS Database DESTINATION ${LIBRARY_OUTPUT_PATH} ) 
