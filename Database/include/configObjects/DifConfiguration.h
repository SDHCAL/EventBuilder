#ifndef _DIFCONFIGURATION_H_
#define _DIFCONFIGURATION_H_

#include "Configuration.h"
#include "Dif.h"
#include <vector>

/**
   Class managing the configuration of DIFs
**/
class DifConfiguration : public Configuration{
 
 private :

 public:
  DifConfiguration();
  ~DifConfiguration();

  void download(int maj, int min) throw (ILCException::Exception);
  /**
     Get the list of unused versions
     @return A vector containing all the unused versions. Once used, you should free the memory using the Version::clearVersionVector(vector<Version*> v) method.
  **/
  static vector<Version*> getUnusedVersions() throw (ILCException::Exception);
};

#endif
