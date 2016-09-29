#################################
# cmake file for building DATABASE
# @author LAGARDE François
#################################

# project name
PROJECT( Analysis )

# project version
SET( ANALYSIS_VERSION_MAJOR 1 )
SET( ANALYSIS_VERSION_MINOR 0 )
SET( ANALYSIS_VERSION_PATCH 0 )
# -------------------------------------------------
# build the analysis library
ADD_DEFINITIONS(-std=c++1y)
file(COPY ./include/ DESTINATION ${PROJECT_INCLUDE_DIR})
if(HaveOCCI)
  FILE (GLOB_RECURSE ANALYSIS_SRC_FILES ./src/*)
  ADD_LIBRARY( Analysis SHARED ${ANALYSIS_SRC_FILES} )
  TARGET_LINK_LIBRARIES( Analysis Reader Database HistoHandler)
  INSTALL(TARGETS Analysis LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib)
else()
  FILE (GLOB_RECURSE ANALYSIS_SRC_FILES ./src/*)
  ADD_LIBRARY( Analysis SHARED ${ANALYSIS_SRC_FILES} )
  TARGET_LINK_LIBRARIES( Analysis Reader HistoHandler)
  INSTALL(TARGETS Analysis LIBRARY  DESTINATION ${CMAKE_SOURCE_DIR}/lib)
endif()

