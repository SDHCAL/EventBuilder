########################################################
# cmake file for building Marlin example Package
# @author Gerald Grenier, IPNL (modifications)
# @author Jan Engels, Desy IT (original)
########################################################





### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )


# add library
AUX_SOURCE_DIRECTORY( ./Geometry library_sources )
ADD_SHARED_LIBRARY( Geometry ${library_sources} )
AUX_SOURCE_DIRECTORY( ./Reader library_sources )
ADD_SHARED_LIBRARY( Reader ${library_sources} )
INSTALL_SHARED_LIBRARY( Geometry DESTINATION lib )
INSTALL_SHARED_LIBRARY( Reader DESTINATION lib )

