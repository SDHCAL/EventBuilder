#include "configObjects/Dif.h"

Diff::Diff(string t) throw (ILCException::Exception):ConfigObject(){
  addParameter("ID","NUMBER");
  addParameter("NAME","VARCHAR2");
  addParameter("DIF_TYPE","VARCHAR2");
  setString("DIF_TYPE",t);
  createParameters("DIF");
  setInt("ENABLED",1);
}

Diff::Diff(int i, string n, string t) throw (ILCException::Exception):ConfigObject(){
  addParameter("ID","NUMBER");
  addParameter("NAME","VARCHAR2");
  addParameter("DIF_TYPE","VARCHAR2");
  setString("DIF_TYPE",t);
  setInt("ID",i);
  setString("NAME",n);
  createParameters("DIF");
  setInt("ENABLED",1);
}

Diff::~Diff(){}

string Diff::toXML() const throw (ILCException::Exception) {
  string s = "<DIF>\n"+ConfigObject::toXML()+"</DIF>";
  return s;
}
