#######################################################
# cmake file for building EventBuilder package
# @author Eté Rémi,Lagarde François IPNL
# Copyright (c) CNRS / IPNL
#######################################################
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-locci" LOCCI_SUPPORT)
if(LOCCI_SUPPORT)
  SET(HaveOCCI 1)
else()
  SET(HaveOCCI 0)
endif()

# configure a header file to pass some of the CMake settings
# to the source code
configure_file ("${PROJECT_SOURCE_DIR}/source/include/Version.h.in" "${PROJECT_BINARY_DIR}/Version.h")
configure_file ("${PROJECT_SOURCE_DIR}/source/include/HaveOCCI.h.in" "${PROJECT_BINARY_DIR}/HaveOCCI.h") 
# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")


# include directories
INCLUDE_DIRECTORIES( include )

# require proper c++
ADD_DEFINITIONS( "-pedantic -Wunused-value -O3 -Wno-long-long -Wreturn-type " )




# -------------------------------------------------
# build the database library
if(HaveOCCI)
FILE (GLOB_RECURSE DATABASE_SRC_FILES src/Database/*)
ADD_SHARED_LIBRARY( Database ${DATABASE_SRC_FILES} )
TARGET_LINK_LIBRARIES( Database)
INSTALL(TARGETS Database LIBRARY DESTINATION lib) 
endif()
# -------------------------------------------------
# build the common library

FILE (GLOB_RECURSE COMMON_SRC_FILES src/Common/*)
ADD_SHARED_LIBRARY( Common ${COMMON_SRC_FILES} )
TARGET_LINK_LIBRARIES( Common )
INSTALL( TARGETS Common LIBRARY DESTINATION lib)

# -------------------------------------------------
# build the analysis library
if(HaveOCCI)
  FILE (GLOB_RECURSE ANALYSIS_SRC_FILES src/Analysis/*)
  ADD_SHARED_LIBRARY( Analysis ${ANALYSIS_SRC_FILES} )
  TARGET_LINK_LIBRARIES( Analysis Common Database )
  INSTALL(TARGETS Analysis LIBRARY DESTINATION lib)
else()
  FILE (GLOB_RECURSE ANALYSIS_SRC_FILES src/Analysis/*)
  ADD_SHARED_LIBRARY( Analysis ${ANALYSIS_SRC_FILES} )
  TARGET_LINK_LIBRARIES( Analysis Common)
  INSTALL(TARGETS Analysis LIBRARY DESTINATION lib)
endif()
# -------------------------------------------------
# build the noise library

FILE (GLOB_RECURSE NOISE_SRC_FILES src/Noise/*)
ADD_SHARED_LIBRARY( Noise ${NOISE_SRC_FILES} )
TARGET_LINK_LIBRARIES( Noise Common )

INSTALL( TARGETS Noise LIBRARY DESTINATION lib)


# -------------------------------------------------
# build the streamout library

FILE (GLOB_RECURSE STREAMOUT_SRC_FILES src/StreamOut/*)

ADD_SHARED_LIBRARY( StreamOut ${STREAMOUT_SRC_FILES} )
#TARGET_LINK_LIBRARIES( StreamOut Common )

INSTALL( TARGETS StreamOut LIBRARY DESTINATION lib
)

# -------------------------------------------------
# build the trivent library

FILE (GLOB_RECURSE TRIVENT_SRC_FILES src/Trivent/*)

ADD_SHARED_LIBRARY( Trivent ${TRIVENT_SRC_FILES} )
TARGET_LINK_LIBRARIES( Trivent Common )

INSTALL(TARGETS Trivent LIBRARY DESTINATION lib)  
