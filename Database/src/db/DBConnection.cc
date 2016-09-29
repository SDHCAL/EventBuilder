#include "db/DBConnection.h"

DBConnection* DBConnection::handle = 0;

DBConnection::DBConnection(){
  user="";
  database="";
  password="";
  env=nullptr;
  conn=nullptr;
}

DBConnection::~DBConnection(){
}

DBConnection* DBConnection::getConnection(){
	if(DBConnection::handle==0){
		DBConnection::handle = new DBConnection();
	}
	return DBConnection::handle;
}

void DBConnection::setUser(string u){
  user = u;
}

void DBConnection::setPassword(string p){
  password = p;
}

void DBConnection::setDatabase(string d){
  database = d;
}

bool DBConnection::isConnected(){
  return (conn!=nullptr);
}

void DBConnection::connect(){
    if(!isConnected()){
      env=Environment::createEnvironment(Environment::DEFAULT);
      conn=env->createConnection(user,password,database); 
    }
  
}

void DBConnection::disconnect(){
    if(isConnected()){
      env->terminateConnection(conn);
      Environment::terminateEnvironment(env);
      env=nullptr;
      conn=nullptr;
    }

 
}

Statement* DBConnection::createStatement(string q){

    if(isConnected()){
      Statement *stmt = conn->createStatement(q);
      return stmt;
    }
    else{
      string m = "DBConnection::createStatement : Not connected to a database.";
      throw ILCException::Exception(m);
    }

  
}

void DBConnection::terminateStatement(Statement* s){
  conn->terminateStatement(s);
}

void DBConnection::commit(){

    conn->commit();

  
}

void DBConnection::rollback(){

    conn->rollback();

  
}

Clob* DBConnection::createClob(){
  Clob* clob = new Clob(conn);
  clob->setEmpty();
  return clob;
}
