#include "configObjects/LdaConfiguration.h"
#include <sstream>

LdaConfiguration::LdaConfiguration():Configuration(){
  objectName="LDA";
}

LdaConfiguration::~LdaConfiguration(){
  clear();
}

void LdaConfiguration::download(int maj, int min) throw (ILCException::Exception){
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

    stmt = db->createStatement("select L.id, C.majorVersionId, C.minorVersionId, L.name, C.firmware, C.hardware from LDA_CONFIG C, LDA L WHERE L.ID=C.LDA_ID AND C.majorVersionId=:max AND C.minorVersionId<=:min order by L.id, C.majorVersionId, C.minorVersionId DESC");
    stmt->setInt(1,maj);
    stmt->setInt(2,min);
    rs = stmt->executeQuery();
    string previousID = "";
    while(rs->next()){
      string id = rs->getString(1);
      if(previousID!=id){
	string n = rs->getString(4);
	string f = rs->getString(5);
	string h = rs->getString(6);
	Lda* l = new Lda(id,n);
	l->setString("FIRMWARE",f);
	l->setString("HARDWARE",h);
	add(l);
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

vector<Version*> LdaConfiguration::getUnusedVersions() throw (ILCException::Exception){
  try{
    return Configuration::getUnusedVersions("lda");
  }
  catch(ILCException::Exception e){
    throw e;
  }
}
