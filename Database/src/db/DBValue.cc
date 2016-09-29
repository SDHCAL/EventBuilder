#include "db/DBValue.h"

string DBValue::getString(){
  throw ILCException::Exception("This value is not a string! (the actual type is "+getType()+")");
}

int DBValue::getInt(){
  throw ILCException::Exception("This value is not an integer (the actual type is "+getType()+")");
}

vector<int> DBValue::getIntVector(){
  throw ILCException::Exception("This value is not a vector of integer (the actual type is "+getType()+")");
}

void DBValue::setString(string s){
  throw ILCException::Exception("This value is not a string (the actual type is "+getType()+")");
}

void DBValue::setInt(int i){
  throw ILCException::Exception("This value is not a integer (the actual type is "+getType()+")");
}

void DBValue::setIntVector(vector<int> v){
  throw ILCException::Exception("This value is not a vector of integer (the actual type is "+getType()+")");
}

void DBValue::setNull(){
  isNull=true;
}

bool DBValue::null(){
  return isNull;
}
