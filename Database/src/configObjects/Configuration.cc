#include "configObjects/Configuration.h"

Configuration::Configuration(){
  modified = false;
}

Configuration::~Configuration(){

}

Version* Configuration::createMajorVersion(string n) const throw (ILCException::Exception){
  Version *newVersion = nullptr;

    DBConnection *db = DBConnection::getConnection();  
    Statement *stmt = db->createStatement("select max(majorVersionId) from "+objectName+"_VERSION");
    ResultSet *rs = stmt->executeQuery();
    if(rs->next()){
      int majorId = rs->getInt(1);
      newVersion = new Version(majorId+1,0,n,false);
      newVersion->create(objectName+"_VERSION");
    }

  
  return newVersion;
}

Version* Configuration::createMinorVersion(int maj, string n) const throw (ILCException::Exception){
  Version *newVersion = nullptr;

    DBConnection *db = DBConnection::getConnection();  
    Statement *stmt = db->createStatement("select max(minorVersionId) from "+objectName+"_VERSION where majorVersionId=:major");
    stmt->setInt(1,maj);
    ResultSet *rs = stmt->executeQuery();
    if(rs->next()){
      int minorId = rs->getInt(1);
      newVersion = new Version(maj,minorId+1,n,false);
      newVersion->create(objectName+"_VERSION");
    }

  
  return newVersion;
}


vector<Version*> Configuration::getVersions(string objectName) throw (ILCException::Exception){
  vector<Version*> v;

    DBConnection *db = DBConnection::getConnection();  
    Statement *stmt = db->createStatement("select majorVersionId, minorVersionId, TO_CHAR(creation,'DD/MM/YYYY HH24:MI:SS'), name, locked from "+objectName+"_VERSION order by majorVersionId, minorVersionId");
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      int majorId = rs->getInt(1);
      int minorId = rs->getInt(2);
      string creation = rs->getString(3);
      string name = rs->getString(4);
      int locked = rs->getInt(5);
      Version *newVersion = new Version(majorId,minorId,creation,name,(locked==1));
      v.push_back(newVersion);
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

  
  return v;
}

Version* Configuration::getVersion(string objectName, int majorId, int minorId) throw (ILCException::Exception){
  Version* v=nullptr;

    DBConnection *db = DBConnection::getConnection();  
    Statement *stmt = db->createStatement("select majorVersionId, minorVersionId, TO_CHAR(creation,'DD/MM/YYYY HH24:MI:SS'), name, locked from "+objectName+"_VERSION WHERE majorVersionId=:m AND minorVersionId=:min");
    stmt->setInt(1,majorId);
    stmt->setInt(2, minorId);
    ResultSet *rs = stmt->executeQuery();
    if(rs->next()){
      int majorId = rs->getInt(1);
      int minorId = rs->getInt(2);
      string creation = rs->getString(3);
      string name = rs->getString(4);
      int locked = rs->getInt(5);
      v = new Version(majorId,minorId,creation,name,(locked==1));
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

  
  return v;
}

void Configuration::clear(){
  for(uint i=0;i<list.size();i++){
    delete(list.at(i));
  }
  list.clear();
}

int Configuration::size(){
  return list.size();
}

void Configuration::add(ConfigObject* l){
  list.push_back(l);
}

vector<ConfigObject*>& Configuration::getVector(){
  return list;
}

void Configuration::setModified(bool b){
  modified = b;
}

bool Configuration::isModified(){
  return modified;
}

bool Configuration::versionExist(string objectName, string n) throw (ILCException::Exception){

     bool result = false;
     DBConnection *db = DBConnection::getConnection();  
     Statement *stmt = db->createStatement("select majorVersionId from "+objectName+"_VERSION WHERE name=:n");
     stmt->setString(1, n);
     ResultSet *rs = stmt->executeQuery();
     if(rs->next()){
       result = true;
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
    return result;

  
}

bool Configuration::minorUpdateNeeded(){
  for(uint i=0;i<list.size();i++){
    if(list[i]->isModified())
      return true;
  }
  return false;
}

vector<Version*> Configuration::getUnusedVersions(string object_name) throw (ILCException::Exception){
  Version* v=nullptr;
  vector<Version*> list;
  
    DBConnection *db = DBConnection::getConnection();  
    Statement *stmt = db->createStatement("select VT.majorVersionId, VT.minorVersionId, TO_CHAR(VT.creation,'DD/MM/YYYY HH24:MI:SS'), VT.name, VT.locked from "+object_name+"_version VT, (select majorVersionId from "+object_name+"_version minus select AV.majorVersionId from "+object_name+"_version AV, state S where AV.majorVersionId=S."+object_name+"majversion) CC WHERE VT.majorVersionId=CC.majorVersionId");
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      int majorId = rs->getInt(1);
      int minorId = rs->getInt(2);
      string creation = rs->getString(3);
      string name = rs->getString(4);
      int locked = rs->getInt(5);
      v = new Version(majorId,minorId,creation,name,(locked==1));
      list.push_back(v);
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  
  
  return list;
}

vector<string> Configuration::toXML(){
  vector<string> v;
  for(vector<ConfigObject*>::iterator it=list.begin();it!=list.end();it++){
    v.push_back((*it)->toXML());
  }
  return v;
}

void Configuration::saveToXML(string fileName){
  vector<string> conf = this->toXML();
  fstream f;
  bool exists=false;
  f.open(fileName.c_str(),ios::in);
  if( f.is_open() ){
    exists=true;
  }
  f.close();

  if(exists){
    f.open(fileName.c_str(),ios::out | ios::app);
    for(vector<string>::iterator it=conf.begin();it!=conf.end();it++){
      f<<(*it)<<endl;
    }
    f.close();
  }
  else{
    f.open(fileName.c_str(),ios::out);
    f<<"<CONFIG>"<<endl;
    for(vector<string>::iterator it=conf.begin();it!=conf.end();it++){
      f<<(*it)<<endl;
    }
    f<<"</CONFIG>"<<endl;
    f.close();
  }
}
