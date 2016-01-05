#include "Database/daq/RunInfo.h"

RunInfo::RunInfo(int t, string d) throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();

  number = -1;
  daqID=-1;
  status = 0;
  type = t;
  startDate="";
  stopDate="";
  comment = "";


    Statement* stmt = db->createStatement("SELECT run_seq.nextval, ID FROM DAQ WHERE NAME=:n");
    stmt->setString(1,d);
    ResultSet* rs = stmt->executeQuery();
    if(rs->next()){
      number = rs->getInt(1);
      daqID = rs->getInt(2);
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

    if(number==-1){
      ostringstream oss;
      oss<<"Can not create run : the DAQ "<<d<<" does not exist!"<<endl;
      throw ILCException::Exception(oss.str());
    }

    stmt = db->createStatement("INSERT INTO RUNINFO (ID, STATUS, RUN_TYPE, START_TIME, STOP_TIME, DAQ_ID, DESCRIPTION) VALUES (:i, :s, :t, :sta, :st, :d, :dd)");
    stmt->setInt(1,number);
    stmt->setInt(2,status);
    stmt->setInt(3,type);
    stmt->setNull(4, OCCIDATE);
    stmt->setNull(5, OCCIDATE);
    stmt->setInt(6,daqID);
    stmt->setString(7,comment);

    stmt->executeUpdate();
    db->terminateStatement(stmt);
    db->commit();

  
}

RunInfo::RunInfo(){
  number=-1;
  status=0;
  type=0;
  startDate="";
  stopDate="";
  daqID=-1;
  comment="";
}

/*void RunInfo::setStatus(int s) throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();

    Statement *stmt = db->createStatement("UPDATE RUNINFO SET STATUS=:s WHERE ID=:n");
    stmt->setInt(1,s);
    stmt->setInt(2,number);
    stmt->executeUpdate();
    db->terminateStatement(stmt);
    status = s;
    db->commit();

  
}*/

/*void RunInfo::start() throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();
  if(startDate.compare("")!=0){
    ostringstream oss;
    oss<<"The run has already been started!";
    throw ILCException::Exception(oss.str());
  }

    Statement *stmt = db->createStatement("UPDATE RUNINFO SET START_TIME=SYSDATE WHERE ID=:n AND START_TIME IS NULL");
    stmt->setInt(1,number);
    stmt->executeUpdate();
    db->terminateStatement(stmt);

    stmt = db->createStatement("SELECT TO_CHAR(START_TIME,'DD/MM/YYYY HH24:MI:SS') FROM RUNINFO WHERE ID=:n");
    stmt->setInt(1,number);
    ResultSet* rs = stmt->executeQuery();
    if(rs->next()){
      startDate = rs->getString(1);
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

    setStatus(2);//running

    db->commit();

  
}*/

/*void RunInfo::stop() throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();

  if(startDate.compare("")==0){
    ostringstream oss;
    oss<<"The run has not been started!";
    throw ILCException::Exception(oss.str());
  }

  if(stopDate.compare("")!=0){
    ostringstream oss;
    oss<<"The run is already stopped!";
    throw ILCException::Exception(oss.str());
  }


    Statement *stmt = db->createStatement("UPDATE RUNINFO SET STOP_TIME=SYSDATE WHERE ID=:n AND STOP_TIME IS NULL");
    stmt->setInt(1,number);
    stmt->executeUpdate();
    db->terminateStatement(stmt);

    stmt = db->createStatement("SELECT TO_CHAR(STOP_TIME,'DD/MM/YYYY HH24:MI:SS') FROM RUNINFO WHERE ID=:n");
    stmt->setInt(1,number);
    ResultSet* rs = stmt->executeQuery();
    if(rs->next()){
      stopDate = rs->getString(1);
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

    setStatus(3);//stopped

    db->commit();

  
}*/

/*void RunInfo::setDescription(string d) throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();

    Statement *stmt = db->createStatement("UPDATE RUNINFO SET DESCRIPTION=:d WHERE ID=:n");
    stmt->setString(1,d);
    stmt->setInt(2,number);
    stmt->executeUpdate();
    db->terminateStatement(stmt);
    comment = d;
    db->commit();

  
}*/

int RunInfo::getRunNumber(){
  return number;
}

int RunInfo::getStatus(){
  return status;
}

int RunInfo::getType(){
  return type;
}

string RunInfo::getStartTime(){
  return startDate;
}

string RunInfo::getStopTime(){
  return stopDate;
}

string RunInfo::getDescription(){
  return comment;
}

Daq* RunInfo::getDaq() throw (ILCException::Exception){
  try{
    string n = Daq::getNameFromId(daqID);
    Daq* d = Daq::getDaq(n);
    return d;
  }
  catch(ILCException::Exception e){
    throw e;
  }
}

RunInfo* RunInfo::getRunInfo(int i) throw (ILCException::Exception){
   DBConnection *db = DBConnection::getConnection();  
   RunInfo* d = NULL;

    Statement* stmt = db->createStatement("SELECT ID, STATUS, RUN_TYPE, TO_CHAR(START_TIME,'DD/MM/YYYY HH24:MI:SS'),  TO_CHAR(STOP_TIME,'DD/MM/YYYY HH24:MI:SS'), DAQ_ID, DESCRIPTION FROM RUNINFO WHERE ID=:n");
    stmt->setInt(1,i);
    ResultSet *rs = stmt->executeQuery();
    if(rs->next()){
      d = new RunInfo();
      d->number = rs->getInt(1);
      d->status = rs->getInt(2);
      d->type = rs->getInt(3);
      d->startDate = rs->getString(4);
      d->stopDate = rs->getString(5);
      d->daqID = rs->getInt(6);
      d->comment = rs->getString(7);
    }
    else{
      ostringstream oss;
      oss<<"There is no run with number "<<i<<"!";
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);      
      throw ILCException::Exception(oss.str());
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

    return d;

  
}

string RunInfo::getStatusDefinition(int val) throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();

    Statement* stmt = db->createStatement("SELECT DESCRIPTION FROM RUN_STATUS WHERE ID=:i");
    stmt->setInt(1,val);
    ResultSet* rs = stmt->executeQuery();
    if(!rs->next()){
      ostringstream oss;
      oss<<"There is no run status "<<val<<"!"<<endl;
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);
      throw ILCException::Exception(oss.str());
    }
    string ret = rs->getString(1);
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
    return ret;

   
}

string RunInfo::getTypeDefinition(int val) throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();

    Statement* stmt = db->createStatement("SELECT DESCRIPTION FROM RUN_TYPE WHERE ID=:i");
    stmt->setInt(1,val);
    ResultSet* rs = stmt->executeQuery();
    if(!rs->next()){
      ostringstream oss;
      oss<<"There is no run type "<<val<<"!"<<endl;
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);
      throw ILCException::Exception(oss.str());
    }
    string ret = rs->getString(1);
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
    return ret;

   
}

vector<RunInfo*> RunInfo::getRunsFromDaqName(string dname) throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();
  vector<RunInfo*> list;

    Statement* stmt = db->createStatement("SELECT R.ID FROM DAQ D, RUNINFO R WHERE D.ID=R.DAQ_ID AND D.NAME=:n ORDER BY R.ID");
    stmt->setString(1,dname);
    ResultSet* rs = stmt->executeQuery();
    while(rs->next()){
      int ret = rs->getInt(1);
      RunInfo* r=RunInfo::getRunInfo(ret);
      list.push_back(r);
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

  
  return list;
}
