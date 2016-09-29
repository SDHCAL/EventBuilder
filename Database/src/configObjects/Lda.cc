#include "configObjects/Lda.h"

Lda::Lda() throw (ILCException::Exception):ConfigObject(){
  addParameter("MAC_ADRESS","VARCHAR2");
  addParameter("NAME","VARCHAR2");
  createParameters("LDA");
}

Lda::Lda(string m, string n) throw (ILCException::Exception):ConfigObject(){
  addParameter("MAC_ADRESS","VARCHAR2");
  addParameter("NAME","VARCHAR2");
  createParameters("LDA");
  setString("MAC_ADRESS",m);
  setString("NAME",n);
}

Lda::~Lda(){}

string Lda::toXML() const throw (ILCException::Exception){
  string s = "<LDA>\n"+ConfigObject::toXML()+"</LDA>";
  return s;
}
