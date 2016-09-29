########################################################
# cmake file for building Marlin
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
########################################################

PROJECT( Marlino )
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 8 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

file(COPY ./include/ DESTINATION ${PROJECT_INCLUDE_DIR})
configure_file ("./include/cmake/Version.h.in" "./include/cmake/Version.h")
file(COPY include/cmake/Version.h DESTINATION ${PROJECT_INCLUDE_DIR})
AUX_SOURCE_DIRECTORY( ./src library_sources )

# need to remove ProcessorLoader.cc and Marlin.cc from library sources
LIST( REMOVE_ITEM library_sources ./src/marlino/Marlino.cc ./src/marlino/ProcessorLoader.cc )

# create library
ADD_LIBRARY( Marlino ${library_sources} )
INSTALL(TARGETS Marlino DESTINATION ${LIBRARY_OUTPUT_PATH} )
TARGET_LINK_LIBRARIES( Marlino lcio ${CMAKE_DL_LIBS} Tinyxml )

# ----- Marlin executable ----------------------------------------------------
ADD_EXECUTABLE( bin_Marlino ./src/Marlino.cc ./src/ProcessorLoader.cc )
# avoid name conflict between Marlin library and Marlin application
SET_TARGET_PROPERTIES( bin_Marlino PROPERTIES OUTPUT_NAME Marlino )
TARGET_LINK_LIBRARIES( bin_Marlino Marlino )
INSTALL( TARGETS bin_Marlino DESTINATION ${CMAKE_BINARY_DIR} )
