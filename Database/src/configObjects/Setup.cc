#include "configObjects/Setup.h"

Setup::Setup(string n){
  name = n;
  locked=false;
}

Setup::Setup(Setup& s){
  name = s.name;
  locked=s.locked;
  for(uint i=0;i<s.states.size();i++){
    states[i]=s.states[i]->clone();
  }
}

Setup::~Setup(){
  for(uint i=0;i<states.size();i++){
    State* s = states[i];
    delete s;
  }
  states.clear();
}

void Setup::addState(State* s){
  states.push_back(s->clone());
}

bool Setup::isLocked(){
  return locked;
}

Setup* Setup::getSetup(string n) throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();  
  Setup* s=nullptr;
  try{
    int id=-1;
    Statement *stmt = db->createStatement("select id, name, locked FROM SETUP WHERE name=:n");
    stmt->setString(1,n);
    ResultSet *rs = stmt->executeQuery();
    if(rs->next()){
      s = new Setup(n);
      s->locked = (rs->getInt(3)==1);
      id = rs->getInt(1);
    }
    else{ 
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);
      ostringstream oss;
      oss<<"Setup "<<n<<" not found!";
      throw ILCException::Exception(oss.str());
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

    if(s!=nullptr){
      // ajout des states
      stmt = db->createStatement("select STA.name FROM STATE STA, SETUP_STATE SS WHERE SS.state_id=STA.id AND SS.setup_id=:i");
      stmt->setInt(1,id);
      rs = stmt->executeQuery();
      while(rs->next()){
	s->states.push_back(State::download(rs->getString(1)));
      }
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);
    }
  }
  
  catch(ILCException::Exception e){
    throw e;
  }
  return s;
}

vector<string> Setup::getSetupNames() throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();  
  vector<string> names;
  try{
    Statement *stmt = db->createStatement("select name FROM SETUP order by name");
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      names.push_back(rs->getString(1));
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
 
  catch(ILCException::Exception e){
    throw e;
  }
  return names;
}

vector<string> Setup::getUnusedSetupNames() throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();  
  vector<string> names;
  try{
    Statement *stmt = db->createStatement("select name from Setup where id not in (select distinct setup_id from run) order by name");
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      names.push_back(rs->getString(1));
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
 
  catch(ILCException::Exception e){
    throw e;
  }
  return names;
}

vector<State*> Setup::getStates(){
  /*
  vector<State*> v;
  for(uint i=0;i<states.size();i++){
    v.push_back(states[i]->clone());
  }
  */
  return states;
}

string Setup::getName(){
  return name;
}

void Setup::lock(){

    DBConnection *db = DBConnection::getConnection();
    Statement *stmt = db->createStatement("UPDATE SETUP SET locked=1 WHERE name=:n");
    stmt->setString(1, name);
    stmt->executeUpdate();
    db->terminateStatement(stmt);
    locked=true;

  
}

bool Setup::isValid(){
  //We check if the same hardware is not used twice in the setup (2 states using the same hardware)

  // Check the LDAs
  set<string> ldaSet;
  for(uint i=0;i<states.size();i++){
    vector<string> ldas = states[i]->getLdaIDs();
    for(uint j=0;j<ldas.size();j++){
      pair<set<string>::iterator,bool> p = ldaSet.insert(ldas[j]);
      if(!p.second){
	ostringstream oss;
	oss<<"The LDA "<<ldas[j]<<" is used several times in the setup "<<name<<"! The setup can not be uploaded to the database!";
	throw ILCException::Exception(oss.str());	
      }
    }
  }

  // Check the DIFs
  set<int> difSet;
  for(uint i=0;i<states.size();i++){
    vector<int> difs = states[i]->getDifIDs();
    for(uint j=0;j<difs.size();j++){
      pair<set<int>::iterator,bool> p = difSet.insert(difs[j]);
      if(!p.second){
	ostringstream oss;
	oss<<"The DIF "<<difs[j]<<" is used several times in the setup "<<name<<"! The setup can not be uploaded to the database!";
	throw ILCException::Exception(oss.str());	
      }
    }
  }

  // The DCC and ASIC are respectively linked to a LDA or a DIF so there is no need to check them if the LDA and DIF are ok. 

  return true;
}

bool Setup::setupExist(){

     bool result = false;
     DBConnection *db = DBConnection::getConnection();  
     Statement *stmt = db->createStatement("select id from SETUP WHERE name=:n");
     stmt->setString(1, name);
     ResultSet *rs = stmt->executeQuery();
     if(rs->next()){
       result = true;
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
    return result;

  
}

Setup* Setup::getSetupFromRun(int runN) throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();
  try{
    string setupName="";
    Statement* s = db->createStatement("SELECT SETUP.name FROM SETUP, RUN WHERE SETUP.ID=RUN.setup_id AND RUN.ID=:n");
    s->setInt(1,runN);
    ResultSet *rs = s->executeQuery();
    if(rs->next()){
      setupName=rs->getString(1);
    }
    s->closeResultSet(rs);	
    db->terminateStatement(s);
  
    if(setupName.compare("")==0)
      return nullptr;

    return getSetup(setupName);
  }
  catch(ILCException::Exception e){
    throw e;
  }
  
}
