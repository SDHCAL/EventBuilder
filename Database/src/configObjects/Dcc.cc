#include "configObjects/Dcc.h"

Dcc::Dcc() throw (ILCException::Exception):ConfigObject(){
  addParameter("ID","NUMBER");
  addParameter("NAME","VARCHAR2");
  addParameter("LDA_ADDRESS","VARCHAR2");
  addParameter("LDA_CHANNEL","NUMBER");
  createParameters("DCC");
}

Dcc::Dcc(int i, string n, string lda_address, int lda_chan) throw (ILCException::Exception):ConfigObject(){
  addParameter("ID","NUMBER");
  addParameter("NAME","VARCHAR2");
  addParameter("LDA_ADDRESS","VARCHAR2");
  addParameter("LDA_CHANNEL","NUMBER");
  createParameters("DCC");
  setInt("ID",i);
  setString("NAME",n);
  setString("LDA_ADDRESS",lda_address);
  setInt("LDA_CHANNEL",lda_chan);
}

Dcc::Dcc(string n, string lda_address, int lda_chan) throw (ILCException::Exception):ConfigObject(){
  addParameter("ID","NUMBER");
  addParameter("NAME","VARCHAR2");
  addParameter("LDA_ADDRESS","VARCHAR2");
  addParameter("LDA_CHANNEL","NUMBER");
  createParameters("DCC");
  setString("NAME",n);
  setString("LDA_ADDRESS",lda_address);
  setInt("LDA_CHANNEL",lda_chan);
}

Dcc::Dcc(string lda_address, int lda_chan) throw (ILCException::Exception):ConfigObject(){
  addParameter("ID","NUMBER");
  addParameter("NAME","VARCHAR2");
  addParameter("LDA_ADDRESS","VARCHAR2");
  addParameter("LDA_CHANNEL","NUMBER");
  createParameters("DCC");
  setString("LDA_ADDRESS",lda_address);
  setInt("LDA_CHANNEL",lda_chan);
}

Dcc::~Dcc(){}

string Dcc::toXML() const throw (ILCException::Exception){
  string s = "<DCC>\n"+ConfigObject::toXML()+"</DCC>";
  return s;
}
