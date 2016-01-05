#ifndef _LDA_H_
#define _LDA_H_

#include <stdio.h>
#include "ConfigObject.h"

/**
   Class representing a LDA
**/
class Lda : public ConfigObject{

 private:

 public:
  Lda() throw (ILCException::Exception);
  Lda(string m, string n) throw (ILCException::Exception);
  ~Lda();

  string getObjectName() const {return "LDA";}
  string toXML() const throw (ILCException::Exception);

};
#endif
