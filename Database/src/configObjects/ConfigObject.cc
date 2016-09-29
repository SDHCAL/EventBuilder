#include "configObjects/ConfigObject.h"

map<string, vector<DBField*> > ConfigObject::dbCache;


ConfigObject::ConfigObject(){
  modified = false;
}

ConfigObject::~ConfigObject(){
  for ( map<string,DBValue*>::iterator it=parameters.begin() ; it != parameters.end(); it++ )
    delete (*it).second;
  for ( map<string,vector<string>*>::iterator it=table_parameter_map.begin() ; it != table_parameter_map.end(); it++ )
    delete (*it).second;
}

void ConfigObject::setModified(bool b){
  modified = b;
}

bool ConfigObject::isModified(){
  return modified;
}

void ConfigObject::addToCache(string n, ResultSet* r){
  if(r==nullptr)
    return;
  vector<DBField*> v;
  while(r->next()){
    v.push_back(new DBField(r->getString(1), r->getString(2), r->getInt(3)));
  }
  dbCache[n]=v;  
}

void ConfigObject::createParameters(string objectName){
  try{
    int go=2;
    string table = objectName;
    while(go>0){
      map<string, vector<DBField*> >::const_iterator it = dbCache.find(table);
      if(it==dbCache.end()){
	DBConnection *db = DBConnection::getConnection();
	Statement *stmt = db->createStatement("select column_name, data_type, data_length from ALL_TAB_COLUMNS where table_name=:name");
	stmt->setString(1,table+"_CONFIG");
	ResultSet* rs = stmt->executeQuery();
	addToCache(table, rs);
	stmt->closeResultSet(rs);
	db->terminateStatement(stmt);
      }
      vector<string>* params = new vector<string>();
      vector<DBField*> vec = dbCache[table];
      for(uint i=0;i<vec.size();i++){
	string name = vec[i]->getName();
	string type = vec[i]->getType();
	int length = vec[i]->getLength();
	
	if(name.compare("MINORVERSIONID")==0 || name.compare("MAJORVERSIONID")==0 || name.compare(table+"_ID")==0){
	  // these are not parameters
	  continue;
	}

	params->push_back(name);
	
	if(type.compare("NUMBER")==0)
	  parameters[name]=new IntVal();
	
	size_t found = type.find("VARCHAR");
	if (found!=string::npos){
	  if(length==1000)//this is a vector of int
	    parameters[name]=new IntVectorVal();
	  else
	    parameters[name]=new StringVal();
	}
      }
      table_parameter_map[table+"_CONFIG"]=params;
      map<string,DBValue*>::const_iterator iter = parameters.find(objectName+"_TYPE");
      if(iter!=parameters.end()){
	table=iter->second->getString();
      }
      else{
	go--;
      }
      go--;
    }
  }
  catch(ILCException::Exception e){cout<<e.getMessage();}
  
}

void ConfigObject::addParameter(string param, string type){
  if(type.compare("NUMBER")==0)
    parameters[param]=new IntVal();
  
  size_t found = type.find("VARCHAR");
  if (found!=string::npos){
    if(type.find("(1000)")!=string::npos)//this is a vector of int
      parameters[param]=new IntVectorVal();
    else
      parameters[param]=new StringVal();
  }
  vector<string>* params = table_parameter_map["BASE"];
  if(params==nullptr)
    params=new vector<string>();
  params->push_back(param);
  table_parameter_map["BASE"]=params;
}

void ConfigObject::setInt(string parameter, int val) throw (ILCException::Exception){
  map<string,DBValue*>::const_iterator iter = parameters.find(parameter);
  if(iter==parameters.end()){
    ostringstream oss;
    oss<<"There is no parameter \""<< parameter <<"\" in the configuration! Available parameters are : "<<endl;
    for ( map<string,DBValue*>::iterator it=parameters.begin() ; it != parameters.end(); it++ )
      oss<<(*it).first<<endl;
    throw ILCException::Exception(oss.str());
  }
  try{
    iter->second->setInt(val);
  }
  catch(ILCException::Exception e){
    throw e;
  }
}

void ConfigObject::setString(string parameter, string val) throw (ILCException::Exception){
  map<string,DBValue*>::const_iterator iter = parameters.find(parameter);
  if(iter==parameters.end()){
    ostringstream oss;
    oss<<"There is no parameter \""<< parameter <<"\" in the configuration! Available parameters are : "<<endl;
    for ( map<string,DBValue*>::iterator it=parameters.begin() ; it != parameters.end(); it++ )
      oss<<(*it).first<<endl;
    throw ILCException::Exception(oss.str());
  }
  try{

      iter->second->setString(val);
  }
  catch(ILCException::Exception e){
    throw e;
  }
}

void ConfigObject::setIntVector(string parameter, vector<int> val) throw (ILCException::Exception){
  map<string,DBValue*>::const_iterator iter = parameters.find(parameter);
  if(iter==parameters.end()){
    ostringstream oss;
    oss<<"There is no parameter \""<< parameter <<"\" in the configuration! Available parameters are : "<<endl;
    for ( map<string,DBValue*>::iterator it=parameters.begin() ; it != parameters.end(); it++ )
      oss<<(*it).first<<endl;
    throw ILCException::Exception(oss.str());
  }
  try{
    iter->second->setIntVector(val);
  }
  catch(ILCException::Exception e){
    throw e;
  }
}

void ConfigObject::setNull(string parameter) throw (ILCException::Exception){
  map<string,DBValue*>::const_iterator iter = parameters.find(parameter);
  if(iter==parameters.end()){
    ostringstream oss;
    oss<<"There is no parameter \""<< parameter <<"\" in the configuration! Available parameters are : "<<endl;
    for ( map<string,DBValue*>::iterator it=parameters.begin() ; it != parameters.end(); it++ )
      oss<<(*it).first<<endl;
    throw ILCException::Exception(oss.str());
  }
  try{
    iter->second->setNull();
  }
  catch(ILCException::Exception e){
    throw e;
  }
}

bool ConfigObject::isNull(string parameter) throw (ILCException::Exception){
  map<string,DBValue*>::const_iterator iter = parameters.find(parameter);
  if(iter==parameters.end()){
    ostringstream oss;
    oss<<"There is no parameter \""<< parameter <<"\" in the configuration! Available parameters are : "<<endl;
    for ( map<string,DBValue*>::iterator it=parameters.begin() ; it != parameters.end(); it++ )
      oss<<(*it).first<<endl;
    throw ILCException::Exception(oss.str());
  }
  try{
    return iter->second->null();
  }
  catch(ILCException::Exception e){
    throw e;
  }
}


int ConfigObject::getInt(string parameter) throw (ILCException::Exception){
  map<string,DBValue*>::const_iterator iter = parameters.find(parameter);
  if(iter==parameters.end()){
    ostringstream oss;
    oss<<"There is no parameter \""<< parameter <<"\" in the configuration! Available parameters are : "<<endl;
    for ( map<string,DBValue*>::iterator it=parameters.begin() ; it != parameters.end(); it++ )
      oss<<(*it).first<<endl;
    throw ILCException::Exception(oss.str());
  }
  try{
    return (iter->second->getInt());
  }
  catch(ILCException::Exception e){
    throw e;
  }
}

string ConfigObject::getString(string parameter) throw (ILCException::Exception){
  map<string,DBValue*>::const_iterator iter = parameters.find(parameter);
  if(iter==parameters.end()){
    ostringstream oss;
    oss<<"There is no parameter \""<< parameter <<"\" in the configuration! Available parameters are : "<<endl;
    for ( map<string,DBValue*>::iterator it=parameters.begin() ; it != parameters.end(); it++ )
      oss<<(*it).first<<endl;
    throw ILCException::Exception(oss.str());
  }
  try{
    return (iter->second->toString());
  }
  catch(ILCException::Exception e){
    throw e;
  }
}

vector<int> ConfigObject::getIntVector(string parameter) throw (ILCException::Exception){
  map<string,DBValue*>::const_iterator iter = parameters.find(parameter);
  if(iter==parameters.end()){
    ostringstream oss;
    oss<<"There is no parameter \""<< parameter <<"\" in the configuration! Available parameters are : "<<endl;
    for ( map<string,DBValue*>::iterator it=parameters.begin() ; it != parameters.end(); it++ )
      oss<<(*it).first<<endl;
    throw ILCException::Exception(oss.str());
  }
  try{
    return (iter->second->getIntVector());
  }
  catch(ILCException::Exception e){
    throw e;
  }
}

string ConfigObject::getType(string parameter) throw (ILCException::Exception){
  map<string,DBValue*>::const_iterator iter = parameters.find(parameter);
  if(iter==parameters.end()){
    ostringstream oss;
    oss<<"There is no parameter \""<< parameter <<"\" in the configuration! Available parameters are : "<<endl;
    for ( map<string,DBValue*>::iterator it=parameters.begin() ; it != parameters.end(); it++ )
      oss<<(*it).first<<endl;
    throw ILCException::Exception(oss.str());
  }
  try{
    return (iter->second->getType());
  }
  catch(ILCException::Exception e){
    throw e;
  }
}



vector<string> ConfigObject::getParametersNames() throw (ILCException::Exception){
  vector<string> list;
  for ( map<string,DBValue*>::iterator it=parameters.begin() ; it != parameters.end(); it++ )
    list.push_back((*it).first);
  return list;
}

vector<string> ConfigObject::getParametersTypes() throw (ILCException::Exception){
  vector<string> list;
  for ( map<string,DBValue*>::iterator it=parameters.begin() ; it != parameters.end(); it++ )
    list.push_back((*it).second->getType());
  return list;
}

vector<string> ConfigObject::getParametersTables() throw (ILCException::Exception){
  vector<string> list;
  for ( map<string, vector<string>*>::iterator it=table_parameter_map.begin() ; it != table_parameter_map.end(); it++ ){
    string tableName = (*it).first;
    if(tableName.substr(0,tableName.find("_CONFIG"))==getObjectName()){
      list.insert(list.begin()+1,tableName);
    }
    else{
      list.push_back(tableName);
    }
  }
  return list;
}

vector<string> ConfigObject::getParametersFromTable(string tableName) throw (ILCException::Exception){
  vector<string> list;
  list = *table_parameter_map[tableName];
  return list;
}


string ConfigObject::cacheToXML(){
  string xml;
  xml.append("<DBPARAMS>\n");
  for ( map<string,vector<DBField*> >::const_iterator it=dbCache.begin() ; it != dbCache.end(); it++ ){
    xml.append("\t<"+it->first+">\n");
    for(unsigned int i=0;i<it->second.size();i++){
      DBField* f = it->second[i];
      xml.append("\t\t<PARAM>");
      xml.append("<NAME>"+f->getName()+"</NAME>");
      xml.append("<TYPE>"+f->getType()+"</TYPE>");
      std::ostringstream oss;
      oss<<f->getLength();
      xml.append("<LENGTH>"+oss.str()+"</LENGTH>");
      xml.append("</PARAM>\n");
    }
    xml.append("\t</"+it->first+">\n");
  }
  xml.append("</DBPARAMS>\n");
  return xml;
}

string ConfigObject::toXML() const{
  string xml;
  for ( map<string,DBValue*>::const_iterator it=parameters.begin() ; it != parameters.end(); it++ )
    xml.append("\t<"+it->first+">"+it->second->toString()+"</"+it->first+">\n");
  return xml;
}

void ConfigObject::terminate(){
  map<string, vector<DBField*> >::const_iterator it;    
  for ( it=dbCache.begin() ; it != dbCache.end(); it++ ){
    vector<DBField*> v = (*it).second;
    for(uint i=0;i<v.size();i++){
      delete(v[i]);
    }
  }
  dbCache.clear();
}

void ConfigObject::setParameters(map<string, string> p){
  map<string,string> param_types;
  vector<string> names = this->getParametersNames();
  vector<string> types = this->getParametersTypes();
  for(uint j=0;j<names.size();j++){
    param_types[names[j]]=types[j];
  }

  for(map<string,string>::const_iterator itr=p.begin();itr!=p.end();itr++){    
    string type = param_types[(*itr).first];
    if(type=="int"){
      int val=0;
      std::istringstream int_val((*itr).second);
      int_val >> val;
      this->setInt((*itr).first,val);
    }
    if(type=="string"){
      this->setString((*itr).first,(*itr).second);
    }
    if(type=="vector<int>"){
      vector<int> values_vec;
      string values_str = (*itr).second;
      char* str=new char [values_str.size()+1];
      strcpy(str,values_str.c_str());
      char* pch=strtok (str," ");
      while (pch != nullptr){
	istringstream ss( pch );
	int val;
	ss >> val;
	values_vec.push_back(val);
	pch = strtok (nullptr, " ");
      }
      delete [] str;
      this->setIntVector((*itr).first,values_vec);
    }
  }
}
