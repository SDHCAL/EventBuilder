#ifndef _DBCONNECTION_H_
#define _DBCONNECTION_H_

#include<occi.h>
#include<stdio.h>
#include<iostream>
#include "Exceptions/Exception.h"

using namespace oracle::occi;
using namespace std;
using namespace ILCException;

/**
   Class used to manage the connection to the database. This is a singleton class, you can not call the constructor. Instead, you must call the DBConnection::getConnection() static method.
**/

class DBConnection{
  
 private:
  string database;
  string user;
  string password;
  Environment *env;
  Connection *conn;

  static DBConnection* handle;
  DBConnection();
  
 public:
  /**
     Get a pointer on the unique instance of the class. Creates the object if needed.
     @return A pointer on the unique DBConnection object
  **/
  static DBConnection* getConnection();
  ~DBConnection();

  /**
     Set the database user name
     @param u The user name
  **/
  void setUser(string u);
 /**
     Set the database user password
     @param p The user password
  **/
  void setPassword(string p);
 /**
     Set the database connection string
     @param d The connection string
  **/
  void setDatabase(string d);

 /**
     Check if we are connected to the database
     @return True if the connection has been performed, false otherwise
  **/
  bool isConnected();
 /**
     Try to connect to the database. Throws a SQLConnection exception in case of problem.
  **/
  void connect();

 /**
     Try to disconnect from the database. Throws a SQLConnection exception in case of problem.
  **/
  void disconnect();

 /**
     Create a statement from the connection. Throws a SQLException or a Exception in case of problem. Once used, you MUST call the terminateStatement method on the statement.
     @param q The query you want to send
     @return A pointer on the newly created statement
  **/
  Statement* createStatement(string q);

  
 /**
    Terminate the given statement.
     @param s A pointer on the statement to close.
  **/
  void terminateStatement(Statement* s);

  /**
     Commit everything since last commit/rollback
  **/
  void commit();

  /**
     Rollback everything since last commit/rollback
  **/
  void rollback();

  /**
     Creates an uninitialized Clob
  **/
  Clob* createClob();

};

#endif
