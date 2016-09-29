#include "daq/Daq.h"

Daq::Daq(string n){
  id = -1;
  status = 0;
  configName = "";
  host = "";
  xml = "";
  comment = "";
  name = n;
  creationDate = "";
}

Daq::~Daq(){
  for(uint i=0;i<jobList.size();i++){
    delete jobList[i];
    jobList.clear();
  }
}

/*void Daq::setConfigName(string s){
  configName=s;
}*/

/*void Daq::setStatus(int s){
  status = s;
}*/

/*void Daq::setHost(string h){
  host = h;
}*/

void Daq::setXML(string fileName) throw (ILCException::Exception){
  string line;
  string text="";
  ifstream myfile(fileName.c_str());
  if (myfile.is_open()){
    while ( myfile.good() ){
      getline (myfile,line);
      text.append(line);
      text.append("\n");
    }
    myfile.close();
  }
  else{
    ostringstream oss;
    oss<<"Unable to open file "<<fileName<<"!"<<endl;
    throw ILCException::Exception(oss.str());
  }
  xml = text;
}

/*void Daq::setComment(string c){
  comment = c;
}*/

string Daq::getConfigName(){
  return configName;
}

int Daq::getStatus(){
  return status;
}

string Daq::getName(){
  return name;
}

string Daq::getTime(){
  return creationDate;
}

string Daq::getHost(){
  return host;
}

string Daq::getXML(){
  return xml;
}

string Daq::getComment(){
  return comment;
}


Daq* Daq::getDaq(string n) throw (ILCException::Exception){
   DBConnection *db = DBConnection::getConnection();  
   Daq* d = nullptr;

    Statement* stmt = db->createStatement("SELECT ID, NAME, STATUS, TO_CHAR(TIME,'DD/MM/YYYY HH24:MI:SS'), CONFIG_NAME, HOST, DESCRIPTION, XML FROM DAQ WHERE NAME=:n");
    stmt->setString(1,n);
    ResultSet *rs = stmt->executeQuery();
    if(rs->next()){
      d = new Daq(rs->getString(2));
      d->id = rs->getInt(1);
      d->status = rs->getInt(3);
      d->creationDate = rs->getString(4);
      d->configName = rs->getString(5);
      d->host = rs->getString(6);
      d->comment = rs->getString(7);
      Clob xmlClob = rs->getClob(8);
      xmlClob.open (OCCI_LOB_READONLY);
      
      Stream *instream = xmlClob.getStream (1,0);
      int size = xmlClob.length();
      char *buffer = new char[size+1];
      memset(buffer, 0, size+1);
      instream->readBuffer (buffer, size);
      string test(buffer);
      d->xml = test;
      delete [] buffer;
      xmlClob.closeStream(instream);
      xmlClob.close();
    }
    else{
      ostringstream oss;
      oss<<"There is no DAQ with name "<<n<<"!";
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);      
      throw ILCException::Exception(oss.str());
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

    stmt = db->createStatement("SELECT ID, JOB_ID, STATUS, DAQ_ID, HOST FROM JOB WHERE DAQ_ID=:i order by host, job_id");
    stmt->setInt(1,d->id);
    rs = stmt->executeQuery();
    while(rs->next()){
      Job* tmpJob = new Job(rs->getString(5),rs->getInt(2));
      tmpJob->setStatus(rs->getInt(3));
      tmpJob->setInternalId(rs->getInt(1));
      d->addJob(tmpJob);
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

    return d;
  
  
}

string Daq::getNameFromId(int i) throw (ILCException::Exception){
   DBConnection *db = DBConnection::getConnection();  

    Statement* stmt = db->createStatement("SELECT NAME FROM DAQ WHERE ID=:i");
    stmt->setInt(1,i);
    ResultSet *rs = stmt->executeQuery();
    if(!rs->next()){
      ostringstream oss;
      oss<<"There is no DAQ with ID "<<i<<"!";
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);      
      throw ILCException::Exception(oss.str());
    }
    else{
      string n = rs->getString(1);
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);
      return n;
    }

 
}


/*void Daq::uploadToDatabase() throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();  
  
    if(xml.length()==0)
      throw ILCException::Exception("Daq::uploadToDatabase() : The XML configuration file should not be empty!");
    Clob xmlClob;
    if(id==-1){ // The object does not exist yet -> insert
      Statement *insertStmt = db->createStatement("INSERT INTO DAQ (ID, NAME, STATUS, TIME, CONFIG_NAME, HOST, DESCRIPTION, XML) VALUES (daq_seq.nextval, :n, :s, SYSDATE, :c, :h, :d, :x)");
      insertStmt->setString(1,name);
      insertStmt->setInt(2,status);
      insertStmt->setString(3,configName);
      insertStmt->setString(4,host);
      insertStmt->setString(5,comment);
      //insert empty clob (it is not possible to directly create a clob)
      Clob* emptyClob = db->createClob();
      insertStmt->setClob(6,*emptyClob);
      insertStmt->executeUpdate();
      db->terminateStatement(insertStmt);

      Statement* stmt = db->createStatement("SELECT ID, TO_CHAR(TIME,'DD/MM/YYYY HH24:MI:SS'), XML FROM DAQ WHERE NAME=:n");
      stmt->setString(1,name);
      ResultSet *rs = stmt->executeQuery();
      if(rs->next()){
	id = rs->getInt(1);
	creationDate = rs->getString(2);
	xmlClob = rs->getClob(3);
      }
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);

      // Conversion from string to clob
      xmlClob.open (OCCI_LOB_READWRITE);
      Stream *outstream = xmlClob.getStream (1,0);
      char* buffer = new char[xml.length()];
      strncpy(buffer, xml.c_str(), xml.length());
      outstream->writeBuffer (buffer, xml.length());
      delete [] buffer;
      outstream->writeLastBuffer ("",0);
      xmlClob.closeStream (outstream);
      xmlClob.close();
      stmt = db->createStatement("UPDATE DAQ SET XML=:xml WHERE ID=:id");
      stmt->setClob(1,xmlClob);
      stmt->setInt(2,id);
      stmt->executeUpdate();
      db->terminateStatement(stmt);
    }
    else{ // the object already exists -> update

      Statement* stmt = db->createStatement("SELECT XML FROM DAQ WHERE ID=:id FOR UPDATE of XML");
      stmt->setInt(1,id);
      ResultSet *rs = stmt->executeQuery();
      if(rs->next()){
	xmlClob = rs->getClob(1);
      }
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);

      Statement *updateStmt = db->createStatement("UPDATE DAQ SET STATUS=:status, CONFIG_NAME=:config, XML=:xml, HOST=:host, DESCRIPTION=:description WHERE ID=:id");
      updateStmt->setInt(1,status);
      updateStmt->setString(2,configName);

      // Conversion from string to clob
      xmlClob.open (OCCI_LOB_READWRITE);
      Stream *outstream = xmlClob.getStream (1,0);
      char* buffer = new char[xml.length()];
      strncpy(buffer, xml.c_str(), xml.length());
      outstream->writeBuffer (buffer, xml.length());
      outstream->writeLastBuffer ("",0);
      xmlClob.closeStream (outstream);
      if(xmlClob.length()>xml.length())
	xmlClob.trim(xml.length());
      xmlClob.close();
      delete [] buffer;
      
      updateStmt->setClob(3,xmlClob);
      updateStmt->setString(4,host);
      updateStmt->setString(5,comment);
      updateStmt->setInt(6,id);
      updateStmt->executeUpdate();
      db->terminateStatement(updateStmt);
    }

    //Update the job list
    for(uint i=0;i<jobList.size();i++){
      Job* j = jobList[i];
      j->upload();
    }


  
  db->commit();
}*/

void Daq::addJob(Job* j){
  j->setDaqId(id);
  jobList.push_back(j);
}

vector<Job*> Daq::getJobs(){
  return jobList;
}

string Daq::getStatusDefinition(int val) throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();
 
    Statement* stmt = db->createStatement("SELECT DESCRIPTION FROM DAQ_STATUS WHERE ID=:i");
    stmt->setInt(1,val);
    ResultSet* rs = stmt->executeQuery();
    if(!rs->next()){
      ostringstream oss;
      oss<<"There is no daq status "<<val<<"!"<<endl;
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);
      throw ILCException::Exception(oss.str());
    }
    string ret = rs->getString(1);
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
    return ret;
  
   
}

vector<string> Daq::getListOfDaqNames() throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();
  vector<string> list;

    Statement* stmt = db->createStatement("SELECT NAME FROM DAQ ORDER BY NAME");
    ResultSet* rs = stmt->executeQuery();
    while(rs->next()){
      string ret = rs->getString(1);
      list.push_back(ret);
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  
   
  return list;
}
