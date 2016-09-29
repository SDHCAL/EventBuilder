#include "configObjects/Asic.h"

Asic::Asic(string t) throw (ILCException::Exception):ConfigObject(){
  addParameter("ID","NUMBER");
  addParameter("ASIC_TYPE","VARCHAR2");
  addParameter("DIF_ID","NUMBER");
  addParameter("NUM","NUMBER");
  setString("ASIC_TYPE",t);
  createParameters("ASIC");
  setInt("ENABLED",1);
}

Asic::Asic(string t, int dif_id, int n) throw (ILCException::Exception):ConfigObject(){
  addParameter("ID","NUMBER");
  addParameter("ASIC_TYPE","VARCHAR2");
  addParameter("DIF_ID","NUMBER");
  addParameter("NUM","NUMBER");
  setString("ASIC_TYPE",t);
  createParameters("ASIC");
  setInt("DIF_ID",dif_id);
  setInt("NUM",n);
  setInt("HEADER",n);
  setInt("ENABLED",1);
}

Asic::~Asic(){}

string Asic::toXML() const throw (ILCException::Exception) {
  string s = "<ASIC>\n"+ConfigObject::toXML()+"</ASIC>";
  return s;
}
