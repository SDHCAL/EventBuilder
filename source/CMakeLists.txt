#######################################################
# cmake file for building EventBuilder package
# @author Eté Rémi,Lagarde François IPNL
# Copyright (c) CNRS / IPNL
#######################################################

# configure a header file to pass some of the CMake settings to the source code
configure_file ("${PROJECT_SOURCE_DIR}/source/include/Version.h.in" "${PROJECT_BINARY_DIR}/Version.h")
# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")


# include directories
INCLUDE_DIRECTORIES( include )

# require proper c++
ADD_DEFINITIONS( "-fPIC -pedantic -Wunused-value -O3 -Wno-long-long -Wreturn-type -Wall" )


# -------------------------------------------------
# build the database library
if(HaveOCCI)
FILE (GLOB_RECURSE DATABASE_SRC_FILES src/Database/*)
ADD_LIBRARY( Database SHARED ${DATABASE_SRC_FILES} )
TARGET_LINK_LIBRARIES( Database)
INSTALL(TARGETS Database LIBRARY DESTINATION lib) 
endif()

# -------------------------------------------------
# build the common library

FILE (GLOB_RECURSE COMMON_SRC_FILES src/Common/*)
ADD_LIBRARY( Common SHARED ${COMMON_SRC_FILES} )
TARGET_LINK_LIBRARIES( Common )
INSTALL( TARGETS Common LIBRARY DESTINATION lib)

# -------------------------------------------------
# build the analysis library

if(HaveOCCI)
  FILE (GLOB_RECURSE ANALYSIS_SRC_FILES src/Analysis/*)
  ADD_LIBRARY( Analysis SHARED ${ANALYSIS_SRC_FILES} )
  TARGET_LINK_LIBRARIES( Analysis Database )
  INSTALL(TARGETS Analysis LIBRARY DESTINATION lib)
else()
  FILE (GLOB_RECURSE ANALYSIS_SRC_FILES src/Analysis/*)
  ADD_LIBRARY( Analysis SHARED ${ANALYSIS_SRC_FILES} )
  TARGET_LINK_LIBRARIES( Analysis Common )
  INSTALL(TARGETS Analysis LIBRARY  DESTINATION lib)
endif()

# -------------------------------------------------
# build the noise library

#FILE (GLOB_RECURSE NOISE_SRC_FILES src/Noise/*)
#ADD_LIBRARY( Noise SHARED ${NOISE_SRC_FILES} )
#TARGET_LINK_LIBRARIES( Noise Common )
#INSTALL( TARGETS Noise LIBRARY DESTINATION lib)


# -------------------------------------------------
# build the streamout library

FILE (GLOB_RECURSE STREAMOUT_SRC_FILES src/StreamOut/*)
ADD_LIBRARY( StreamOut SHARED ${STREAMOUT_SRC_FILES} )
TARGET_LINK_LIBRARIES( StreamOut Common )
INSTALL( TARGETS StreamOut LIBRARY  DESTINATION lib)

# -------------------------------------------------
# build the trivent library

FILE (GLOB_RECURSE TRIVENT_SRC_FILES src/Trivent/*)
ADD_LIBRARY( Trivent SHARED ${TRIVENT_SRC_FILES} )
TARGET_LINK_LIBRARIES( Trivent Common )
INSTALL(TARGETS Trivent LIBRARY DESTINATION lib)  

# ----- ElogFiller executable ----------------------------------------------------
ADD_EXECUTABLE( ElogFiller ${PROJECT_SOURCE_DIR}/source/src/ElogFiller/ElogFiller.cc)
INSTALL( TARGETS ElogFiller DESTINATION bin )
 
# -------------------------------------------------
# build the VirerMarlin library

#FILE (GLOB_RECURSE VIRE_SRC_FILES src/VirerMarlin/*)
#ADD_LIBRARY( VirerMarlin SHARED ${VIRE_SRC_FILES} )
#TARGET_LINK_LIBRARIES( VirerMarlin)
#INSTALL(TARGETS VirerMarlin LIBRARY DESTINATION lib)


# ----- Marlin executable ----------------------------------------------------
#ADD_EXECUTABLE( bin_Marlin 
#${PROJECT_SOURCE_DIR}/source/src/VirerMarlin/Marlin.cc 
#${PROJECT_SOURCE_DIR}/source/src/VirerMarlin/ProcessorLoader.cc 
#)
# avoid name conflict between Marlin library and Marlin application
#SET_TARGET_PROPERTIES( bin_Marlin PROPERTIES OUTPUT_NAME Marlin )
#INSTALL( TARGETS bin_Marlin DESTINATION bin )
# ----------------------------------------------------------------------------


# ----- Simulation executable ----------------------------------------------------
ADD_EXECUTABLE( Simulation ${PROJECT_SOURCE_DIR}/source/src/Simulation/Simulation.cc )
INSTALL( TARGETS Simulation DESTINATION bin )


