#ifndef _ASIC_H_
#define _ASIC_H_

#include <stdio.h>
#include "ConfigObject.h"

/**
   Class representing a ASIC
**/
class Asic : public ConfigObject{

 private:

 public:
  Asic(string t) throw (ILCException::Exception);
  Asic(string t, int dif_id, int n) throw (ILCException::Exception);
  ~Asic();

  string getObjectName() const {return "ASIC";}
  string toXML() const throw (ILCException::Exception);

};
#endif
