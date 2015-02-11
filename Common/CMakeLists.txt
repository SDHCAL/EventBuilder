########################################################
# cmake file for building Marlin example Package
# @author Gerald Grenier, IPNL (modifications)
# @author Jan Engels, Desy IT (original)
########################################################





### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long -std=c++11" )

# include directories
INCLUDE_DIRECTORIES( ./Geometry/include )
INCLUDE_DIRECTORIES( ./Reader/include )
INCLUDE_DIRECTORIES( ./TimeHit/include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# add library
AUX_SOURCE_DIRECTORY( ./Geometry/src library_sources )
ADD_SHARED_LIBRARY( Geometry ${library_sources} )
AUX_SOURCE_DIRECTORY( ./Reader/src library_sources )
ADD_SHARED_LIBRARY( Reader ${library_sources} )
AUX_SOURCE_DIRECTORY( ./TimeHit/src library_sources )
ADD_SHARED_LIBRARY( TimeHit ${library_sources} )
INSTALL_SHARED_LIBRARY( Geometry DESTINATION lib )
INSTALL_SHARED_LIBRARY( Reader DESTINATION lib )
INSTALL_SHARED_LIBRARY( TimeHit DESTINATION lib )


