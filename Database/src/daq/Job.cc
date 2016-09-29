#include "daq/Job.h"

Job::Job(string h, int jid){
  id = -1;
  status = 0;
  host = h;
  jobID=jid;
  daqID=-1;
}

void Job::setStatus(int s){
  status = s;
}

void Job::setInternalId(int i){
  id = i;
}

void Job::setDaqId(int d){
  daqID = d;
}

int Job::getStatus(){
  return status;
}

string Job::getHost(){
  return host;
}

int Job::getID(){
  return jobID;
}

/*void Job::upload() throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();

    if(daqID==-1){
      ostringstream oss;
      oss<<"This job is not associated to a DAQ!";
      throw ILCException::Exception(oss.str());
    }
    if(id==-1){ // The object does not exist yet -> insert
      int newID=-1;
      Statement* idStmt = db->createStatement("SELECT job_seq.nextval FROM DUAL");
      idStmt->executeQuery();
      ResultSet *rs = idStmt->executeQuery();
      if(rs->next()){
	newID = rs->getInt(1);
      }
      idStmt->closeResultSet(rs);
      db->terminateStatement(idStmt);

      Statement *insertStmt = db->createStatement("INSERT INTO JOB (ID, JOB_ID, STATUS, DAQ_ID, HOST) VALUES (:i, :j, :s, :d, :h)");
      insertStmt->setInt(1,newID);
      insertStmt->setInt(2,jobID);
      insertStmt->setInt(3,status);
      insertStmt->setInt(4,daqID);
      insertStmt->setString(5,host);
      insertStmt->executeUpdate();
      db->terminateStatement(insertStmt);

      id = newID;
    }
    else{ // the object already exists -> update

      Statement *updateStmt = db->createStatement("UPDATE JOB SET STATUS=:status, JOB_ID=:j, DAQ_ID=:d, HOST=:h WHERE ID=:id");
      updateStmt->setInt(1,status);
      updateStmt->setInt(2,jobID);
      updateStmt->setInt(3,daqID);
      updateStmt->setString(4,host);
      updateStmt->setInt(5,id);
      updateStmt->executeUpdate();
      db->terminateStatement(updateStmt);
    }

 
  
  //No commit : it will be done in the Daq::uploadToDatabase method!;
}*/

string Job::getStatusDefinition(int val) throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();

    Statement* stmt = db->createStatement("SELECT DESCRIPTION FROM JOB_STATUS WHERE ID=:i");
    stmt->setInt(1,val);
    ResultSet* rs = stmt->executeQuery();
    if(!rs->next()){
      ostringstream oss;
      oss<<"There is no job status "<<val<<"!"<<endl;
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);
      throw ILCException::Exception(oss.str());
    }
    string ret = rs->getString(1);
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
    return ret;

   
}
