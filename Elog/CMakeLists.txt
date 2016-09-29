#################################
# cmake file for building ElogFiller
# @author LAGARDE François
#################################

# project name
PROJECT( Elog )

# project version
SET( ELOG_VERSION_MAJOR 1 )
SET( ELOG_VERSION_MINOR 0 )
SET( ELOG_VERSION_PATCH 0 )
ADD_DEFINITIONS(-std=c++1y)
add_executable(ElogFiller ./src/ElogFiller.cc)
TARGET_LINK_LIBRARIES( ElogFiller Tinyxml )
install(TARGETS ElogFiller DESTINATION ${EXECUTABLE_OUTPUT_PATH})
