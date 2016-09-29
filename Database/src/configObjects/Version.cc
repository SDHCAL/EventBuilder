#include "configObjects/Version.h"
#include <sstream>

Version::Version(int maj, int min, string n, bool l){
  majorId = maj;
  minorId = min;
  name = n;
  creationDate = "";
  locked = l;  
}

Version::Version(int maj, int min, string d, string n, bool l){
  majorId = maj;
  minorId = min;
  name = n;
  creationDate = d;
  locked = l;
}

Version::Version(int maj, int min){
  majorId = maj;
  minorId = min;
  name = "";
  creationDate = "";
  locked = false;
}

Version::~Version(){}

int Version::getMajorId(){
  return majorId;
}

int Version::getMinorId(){
  return minorId;
}

string Version::getName(){
  return name;
}

string Version::getDate(){
  return creationDate;
}

bool Version::isLocked(){
  return locked;
}

void Version::create(string table){

    DBConnection *db = DBConnection::getConnection();
    Statement *stmt = db->createStatement("INSERT INTO "+table+" (majorVersionId, minorVersionId, creation, name, locked) VALUES (:major, :minor, SYSDATE, :name, :locked)");
    stmt->setInt(1, majorId);
    stmt->setInt(2, minorId);
    stmt->setString(3, name);
    stmt->setInt(4, locked);
    stmt->executeUpdate();
    db->terminateStatement(stmt);
    stmt = db->createStatement("select TO_CHAR(creation,'DD/MM/YYYY HH24:MI:SS') from "+table+" WHERE majorVersionId=:majorId AND minorVersionId=:minorId");
    stmt->setInt(1,majorId);
    stmt->setInt(2,minorId);
    ResultSet *rs = stmt->executeQuery();
    if(rs->next()){
      creationDate = rs->getString(1);
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

 
}

void Version::lock(string table){
  
    DBConnection *db = DBConnection::getConnection();
    Statement *stmt = db->createStatement("UPDATE "+table+" SET locked=1 WHERE majorVersionId=:majorId AND minorVersionId=:minorId");
    stmt->setInt(1, majorId);
    stmt->setInt(2, minorId);
    stmt->executeUpdate();
    db->terminateStatement(stmt);
    locked=true;

 
}

void Version::clearVersionVector(vector<Version*> v){
  for(uint i=0;i<v.size();i++){
    delete(v.at(i));
  }
  v.clear();
}

Version* Version::clone(){
  return new Version(majorId, minorId, creationDate, name, locked);
}

string Version::incrementName(string n){
  size_t found;

  //search for the last '_' caracter
  found=n.rfind("_");
  if (found!=string::npos){
    string post = n.substr(found+1, n.length());
    //check if post is an int
    int i;
    istringstream myStream(post);
    if (myStream>>i){
      //post was starting with an int
      //now check that everything was an int...
      ostringstream oss;
      oss<<i;
      if(post.compare(oss.str())==0){
	i++;
	ostringstream oss2;
	oss2<<i;
	return n.substr(0,found+1)+oss2.str();
      }
    }
  }
  return n+"_1";
}

bool Version::isMostRecent(string table){
  bool result = false;

    DBConnection *db = DBConnection::getConnection();
    Statement *stmt = db->createStatement("SELECT * FROM "+table+" WHERE majorVersionId=:majorId AND minorVersionId>:minorId");
    stmt->setInt(1, majorId);
    stmt->setInt(2, minorId);
    ResultSet *rs = stmt->executeQuery();
    if(rs->next()){
      result=false;
    }
    else{
      cout<<"No result found for "<<majorId<<"."<<minorId<<endl;
      result=true;
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
    return result;

  
}
