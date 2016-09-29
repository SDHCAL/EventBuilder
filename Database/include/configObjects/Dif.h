#ifndef _DIF_H_
#define _DIF_H_

#include <stdio.h>
#include "ConfigObject.h"

/**
   Class representing a DIF
**/
class Diff : public ConfigObject{

 private:

 public:
  /**
     Constructor
     @param t The type of the DIF (ex: 'DHCALDIF')
  **/
  Diff(string t) throw (ILCException::Exception);
  /**
     Constructor
     @param id The DIF Num
     @param n The name of the DIF
     @param t The type of the DIF (ex: 'DHCALDIF')
  **/
  Diff(int id, string n, string t) throw (ILCException::Exception);
  ~Diff();

  string getObjectName() const {return "DIF";}
  string toXML() const throw (ILCException::Exception);

};
#endif
