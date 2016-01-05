#ifndef _CONFIGOBJECT_H_
#define _CONFIGOBJECT_H_

#include<stdio.h>
#include<iostream>
#include<cstring>
#include<sstream>
#include <map>
#include "DBField.h"
#include "../db/DBConnection.h"
#include "../db/DBValue.h"
#include "../db/IntVal.h"
#include "../db/IntVectorVal.h"
#include "../db/StringVal.h"

using namespace std;

/**
   Abstract super class for basic objects (LDA, DIF, ASIC, ...)
**/
class ConfigObject{

 private:
  bool modified;
  map<string, DBValue*> parameters;
  map<string, vector<string>* > table_parameter_map;
  
  static map<string, vector<DBField*> > dbCache;//used to avoid a query per object...
  
  void addToCache(string n, ResultSet* r);

 protected:
  void createParameters(string objectName);
  void addParameter(string param, string type);

 public:
  /**
     Returns a XML description of the object
     @return A string containing the XML description
  **/
  virtual string toXML() const;
  /**
     Returns a XML description of the tables descriptions in the cache
     @return A string containing the XML description
  **/
  static string cacheToXML();


  ConfigObject();
  virtual ~ConfigObject();

  /**
     Get the name of the actual object
   **/
  virtual string getObjectName() const=0;

  /**
     Flag the object as beeing modified (or not modified). Only modified objects will be uploaded in a minor version.
     @param b True if the object was modified, false otherwise.
  **/
  void setModified(bool b);
  /**
     Get the modified status of the object
     @return True if the object was said to be modified, false otherwise.
  **/
  bool isModified();

  /**
     Change a int parameter of the object
     @param parameter The name of the parameter.
     @param val The value to set to the parameter
  **/
  void setInt(string parameter, int val) throw (ILCException::Exception);
  /**
     Change a string parameter of the object
     @param parameter The name of the parameter.
     @param val The value to set to the parameter
  **/
  void setString(string parameter, string val) throw (ILCException::Exception);
  /**
     Change a vector<int> parameter of the object
     @param parameter The name of the parameter.
     @param val The value to set to the parameter
  **/
  void setIntVector(string parameter, vector<int> val) throw (ILCException::Exception);
  /**
     Set a parameter of the object to NULL value
     @param parameter The name of the parameter.
  **/
  void setNull(string parameter) throw (ILCException::Exception);

  /**
     Check if the given parameter is null
     @param parameter The name of the parameter.
     @return True if the parameter is null
  **/
  bool isNull(string parameter) throw (ILCException::Exception);
  /**
     Get an int parameter value
     @param parameter The name of the parameter.
     @return The value of the paramter
  **/
  int getInt(string parameter) throw (ILCException::Exception);
  /**
     Get a string parameter value
     @param parameter The name of the parameter.
     @return The value of the paramter
  **/
  string getString(string parameter) throw (ILCException::Exception);
  /**
     Get an vector<int> parameter value
     @param parameter The name of the parameter.
     @return The value of the paramter
  **/
  vector<int> getIntVector(string parameter) throw (ILCException::Exception);
  /**
     Get The type of the given paramter
     @param parameter The name of the parameter.
     @return A string describing the type ('int', 'string', 'vector<int>')
  **/
  string getType(string parameter) throw (ILCException::Exception);

  /**
     Get the list of the parameters names
     @return A vector containing all the parameters names
  **/
  vector<string> getParametersNames() throw (ILCException::Exception);
  /**
     Get the list of the parameters types
     @return A vector containing all the parameters types (the order is the same than in the getParametersNames method)
  **/
  vector<string> getParametersTypes() throw (ILCException::Exception);
  /**
     Get the list of the tables where the parameters are stored
     @return A vector containing all the tables names
  **/
  vector<string> getParametersTables() throw (ILCException::Exception);
  /**
     Get the list of the parameters in the given table
     @return A vector containing all the parameters names
  **/
  vector<string> getParametersFromTable(string tableName) throw (ILCException::Exception);

  static void terminate();

  /**
     Set all parameters contained in the map to the object
     @param p A map containing <parameter_name, parameter_value>
  **/
  void setParameters(map<string, string> p);

};
#endif
