#################################
# cmake file for building DATABASE
# @author LAGARDE François
#################################

# project name
PROJECT( Common )

# project version
SET( COMMON_VERSION_MAJOR 1 )
SET( COMMON_VERSION_MINOR 0 )
SET( COMMON_VERSION_PATCH 0 )

file(COPY ./include/ DESTINATION ${PROJECT_INCLUDE_DIR})
