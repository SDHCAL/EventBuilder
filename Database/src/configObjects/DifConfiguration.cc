#include "configObjects/DifConfiguration.h"
#include <sstream>

DifConfiguration::DifConfiguration():Configuration(){
  objectName="DIF";
}

DifConfiguration::~DifConfiguration(){
  clear();
}

void DifConfiguration::download(int maj, int min) throw (ILCException::Exception){
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

    // Get the different DIF types in the major version
    string query="SELECT DISTINCT DIF_TYPE FROM DIF D, DIF_CONFIG DC WHERE D.ID=DC.DIF_ID AND DC.majorVersionId=:maj AND DC.minorVersionId=0";
    stmt = db->createStatement(query);
    stmt->setInt(1,maj);
    rs = stmt->executeQuery();
    vector<string> dif_types;
    while(rs->next()){
      dif_types.push_back(rs->getString(1));
    }

    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
    
    // we will create 1 query per DIF Type
    for(vector<string>::const_iterator current_type=dif_types.begin();current_type!=dif_types.end();current_type++)
    {
      //Retrieve the parameters of a DIF
      Diff dif(*current_type);
      vector<string> params = dif.getParametersNames();
      vector<string> types = dif.getParametersTypes();
      
      // Create the query from the parameters
      query = "select O.id, C.majorVersionId, C.minorVersionId";
      for(vector<string>::const_iterator itr=params.begin();itr!=params.end();itr++){
	query=query+","+(*itr);
      }
      query+=" FROM DIF_CONFIG C, DIF O, "+(*current_type)+"_CONFIG DC WHERE O.ID=C.DIF_ID AND C.DIF_ID=DC."+(*current_type)+"_ID AND C.MAJORVERSIONID=DC.MAJORVERSIONID AND C.MINORVERSIONID=DC.MINORVERSIONID AND C.majorVersionId=:max AND C.minorVersionId<=:min order by O.id, C.majorVersionId, C.minorVersionId DESC";
      stmt = db->createStatement(query);
      stmt->setInt(1,maj);
      stmt->setInt(2,min);
      // We have many small lines so we prefetch a large number of them to avoid the creation of many connections (takes too much time)
      stmt->setPrefetchMemorySize(0);
      stmt->setPrefetchRowCount(1000);
      rs = stmt->executeQuery();
      int previousID = -1;
      while(rs->next()){
	int id = rs->getInt(1);
	if(previousID!=id){
	  Diff* d = new Diff(*current_type);
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
  }
  
  catch(ILCException::Exception e){
    throw e;
  }
}

vector<Version*> DifConfiguration::getUnusedVersions() throw (ILCException::Exception){
  try{
    return Configuration::getUnusedVersions("dif");
  }
  catch(ILCException::Exception e){
    throw e;
  }
}
