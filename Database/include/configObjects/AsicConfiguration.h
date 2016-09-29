#ifndef _ASICCONFIGURATION_H_
#define _ASICCONFIGURATION_H_

#include "Configuration.h"
#include "Asic.h"
#include <vector>

/**
   Class managing the configuration of ASICs
**/
class AsicConfiguration : public Configuration{
 
 private :

 public:
  AsicConfiguration();
  ~AsicConfiguration();

  void download(int maj, int min) throw (ILCException::Exception);
  /**
     Get the list of unused versions
     @return A vector containing all the unused versions. Once used, you should free the memory using the Version::clearVersionVector(vector<Version*> v) method.
  **/
  static vector<Version*> getUnusedVersions() throw (ILCException::Exception);
};

#endif
