########################################################
# cmake file for building Marlin example Package
# @author Gerald Grenier, IPNL (modifications)
# @author Jan Engels, Desy IT (original)
########################################################





### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic -std=c++0x" )
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ../include )

# add library
AUX_SOURCE_DIRECTORY( . library_sources )
ADD_SHARED_LIBRARY( Analysis ${library_sources} )
INSTALL_SHARED_LIBRARY( Analysis DESTINATION lib )


