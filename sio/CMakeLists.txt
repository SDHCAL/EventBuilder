#################################
# cmake file for building SIO
# @author Frank Gaede, DESY
# @author Jan Engels, DESY
#################################

# project name
PROJECT( SIO )

# project version
SET( SIO_VERSION_MAJOR ${LCIO_VERSION_MAJOR} )
SET( SIO_VERSION_MINOR ${LCIO_VERSION_MINOR} )
SET( SIO_VERSION_PATCH ${LCIO_VERSION_PATCH} )

# library *nix style versioning
SET( SIO_SOVERSION "${SIO_VERSION_MAJOR}.${SIO_VERSION_MINOR}" )
SET( SIO_VERSION "${SIO_SOVERSION}.${SIO_VERSION_PATCH}" )

FIND_PACKAGE( ZLIB REQUIRED )

INCLUDE( TestBigEndian )
TEST_BIG_ENDIAN( BIG_ENDIAN )

IF( NOT BIG_ENDIAN )
    ADD_DEFINITIONS( "-D__LITTLE_ENDIAN__" )
ENDIF()

file(COPY ${SIO_SOURCE_DIR}/include/ DESTINATION ${PROJECT_INCLUDE_DIR})

INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})

SET (SIO_SRCS
    ${SIO_SOURCE_DIR}/src/SIO_block.cc
    ${SIO_SOURCE_DIR}/src/SIO_blockManager.cc
    ${SIO_SOURCE_DIR}/src/SIO_functions.cc
    ${SIO_SOURCE_DIR}/src/SIO_record.cc
    ${SIO_SOURCE_DIR}/src/SIO_recordManager.cc
    ${SIO_SOURCE_DIR}/src/SIO_stream.cc
    ${SIO_SOURCE_DIR}/src/SIO_streamManager.cc
)

ADD_LIBRARY( sio STATIC ${SIO_SRCS} )
TARGET_LINK_LIBRARIES( sio ${ZLIB_LIBRARIES} )
INSTALL( TARGETS sio DESTINATION ${LIBRARY_OUTPUT_PATH} )
