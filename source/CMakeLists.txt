#######################################################
# cmake file for building EventBuilder package
# @author Eté Rémi, IPNL
# Copyright (c) CNRS / IPNL
#######################################################

# include directories
INCLUDE_DIRECTORIES( include )

# require proper c++
ADD_DEFINITIONS( "-pedantic -Wunused-value -O3" )
ADD_DEFINITIONS("-Wno-long-long -Wreturn-type")

# -------------------------------------------------
# build the common library

FILE (
	GLOB_RECURSE
	COMMON_SRC_FILES
	src/Common/*
)

ADD_SHARED_LIBRARY( Common ${COMMON_SRC_FILES} )

INSTALL(
	TARGETS Common
	LIBRARY DESTINATION lib
)

# -------------------------------------------------
# build the analysis library

FILE (
	GLOB_RECURSE
	ANALYSIS_SRC_FILES
	src/Analysis/*
)

ADD_SHARED_LIBRARY( Analysis ${ANALYSIS_SRC_FILES} )
TARGET_LINK_LIBRARIES( Analysis Common )

INSTALL(
	TARGETS Analysis
	LIBRARY DESTINATION lib
)

# -------------------------------------------------
# build the noise library

FILE (
	GLOB_RECURSE
	NOISE_SRC_FILES
	src/Noise/*
)

ADD_SHARED_LIBRARY( Noise ${NOISE_SRC_FILES} )
TARGET_LINK_LIBRARIES( Noise Common )

INSTALL(
	TARGETS Noise
	LIBRARY DESTINATION lib
)


# -------------------------------------------------
# build the streamout library

FILE (
	GLOB_RECURSE
	STREAMOUT_SRC_FILES
	src/StreamOut/*
)

ADD_SHARED_LIBRARY( StreamOut ${STREAMOUT_SRC_FILES} )
#TARGET_LINK_LIBRARIES( StreamOut Common )

INSTALL(
	TARGETS StreamOut
	LIBRARY DESTINATION lib
)

# -------------------------------------------------
# build the trivent library

FILE (
	GLOB_RECURSE
	TRIVENT_SRC_FILES
	src/Trivent/*
)

ADD_SHARED_LIBRARY( Trivent ${TRIVENT_SRC_FILES} )
TARGET_LINK_LIBRARIES( Trivent Common )

INSTALL(
	TARGETS Trivent
	LIBRARY DESTINATION lib
)  
