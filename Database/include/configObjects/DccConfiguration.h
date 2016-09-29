#ifndef _DCCCONFIGURATION_H_
#define _DCCCONFIGURATION_H_

#include "Configuration.h"
#include "Dcc.h"
#include <vector>

/**
   Class managing the configuration of DCCs
**/
class DccConfiguration : public Configuration{
 
 private :

 public:
  DccConfiguration();
  ~DccConfiguration();

  void download(int maj, int min) throw (ILCException::Exception);

  //void addFromXML(string fileName) throw (ILCException::Exception);
  /**
     Get the list of unused versions
     @return A vector containing all the unused versions. Once used, you should free the memory using the Version::clearVersionVector(vector<Version*> v) method.
  **/
  static vector<Version*> getUnusedVersions() throw (ILCException::Exception);
};

#endif
