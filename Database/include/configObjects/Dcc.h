#ifndef _DCC_H_
#define _DCC_H_

#include <stdio.h>
#include "ConfigObject.h"

/**
   Class representing a DCC
**/
class Dcc : public ConfigObject{

 private:

 public:
  Dcc() throw (ILCException::Exception);
  Dcc(int id, string n, string lda_address, int lda_chan) throw (ILCException::Exception);
  Dcc(string n, string lda_address, int lda_chan) throw (ILCException::Exception);
  Dcc(string lda_address, int lda_chan) throw (ILCException::Exception);
  ~Dcc();

  string getObjectName() const {return "DCC";}
  string toXML() const throw (ILCException::Exception);

};
#endif
