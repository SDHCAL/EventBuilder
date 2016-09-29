#include "configObjects/DccConfiguration.h"
#include <sstream>

DccConfiguration::DccConfiguration():Configuration(){
  objectName="DCC";
}

DccConfiguration::~DccConfiguration(){
  clear();
}

void DccConfiguration::download(int maj, int min) throw (ILCException::Exception){
  clear();
  DBConnection *db = DBConnection::getConnection();  
  try{
    //Check the Version does exist
    Statement *stmt = db->createStatement("select majorVersionId from "+objectName+"_VERSION WHERE majorVersionId=:max AND minorVersionId=:min");
    stmt->setInt(1,maj);
    stmt->setInt(2,min);
    ResultSet *rs = stmt->executeQuery();
    if(!rs->next()){
      //no such version
      stmt->closeResultSet(rs); 
      db->terminateStatement(stmt);
      ostringstream oss;
      oss<<"Version "<<maj<<"."<<min<<" does not exist int table "<<objectName<<"_VERSION!";
      throw ILCException::Exception(oss.str());
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

    //Retrieve the parameters of a DCC
    Dcc dcc;
    vector<string> params = dcc.getParametersNames();
    vector<string> types = dcc.getParametersTypes();

    // Create the query from the parameters
    string query = "select O.id, C.majorVersionId, C.minorVersionId";
    for(vector<string>::const_iterator itr=params.begin();itr!=params.end();itr++){
      query=query+","+(*itr);
    }
    query+=" FROM DCC_CONFIG C, DCC O WHERE O.ID=C.DCC_ID AND C.majorVersionId=:max AND C.minorVersionId<=:min order by O.id, C.majorVersionId, C.minorVersionId DESC";
    stmt = db->createStatement(query);
    stmt->setInt(1,maj);
    stmt->setInt(2,min);
    rs = stmt->executeQuery();
    int previousID = -1;
    while(rs->next()){
      int id = rs->getInt(1);
      if(previousID!=id){
	Dcc* d = new Dcc();
	int cpt = 4;
	for(uint i=0;i<params.size();i++){
	  if(rs->isNull(cpt)){
	    d->setNull(params[i]);
	  }
	  else{
	    if(types[i]=="int")
	      d->setInt(params[i],rs->getInt(cpt));
	    if(types[i]=="string")
	      d->setString(params[i],rs->getString(cpt));
	    if(types[i]=="vector<int>")
	      d->setString(params[i],rs->getString(cpt));
	  }
	  cpt++;
	}
	add(d);
      }
      previousID=id;
    }

    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
  
  catch(ILCException::Exception e){
    throw e;
  }
}

vector<Version*> DccConfiguration::getUnusedVersions() throw (ILCException::Exception){
  try{
    return Configuration::getUnusedVersions("dcc");
  }
  catch(ILCException::Exception e){
    throw e;
  }
}

