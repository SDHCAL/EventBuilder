#include "configObjects/State.h"

State::State(){
  name = "UNKNOWN";
  lda_version = nullptr;
  dcc_version = nullptr;
  dif_version = nullptr;
  asic_version = nullptr;
  downloadedLdaConf=nullptr;
  downloadedDccConf=nullptr;
  downloadedDifConf=nullptr;
  downloadedAsicConf=nullptr;
  locked=false;
  comment="";
  nodb = false;
}

State::State(string n){
  name = n;
  lda_version = nullptr;
  dcc_version = nullptr;
  dif_version = nullptr;
  asic_version = nullptr;
  downloadedLdaConf=nullptr;
  downloadedDccConf=nullptr;
  downloadedDifConf=nullptr;
  downloadedAsicConf=nullptr;
  locked=false;
  comment="";
  nodb = false;
}

State::State(const State& s){
  name = s.name;
  if(s.lda_version!=nullptr)
    lda_version = s.lda_version->clone();
  else
    lda_version=nullptr;
  if(s.dcc_version!=nullptr)
    dcc_version = s.dcc_version->clone();
  else
    dcc_version=nullptr;
  if(s.dif_version!=nullptr)
    dif_version = s.dif_version->clone();
  else
    dif_version=nullptr;
  if(s.asic_version!=nullptr)
    asic_version = s.asic_version->clone();
  else
    asic_version=nullptr;
  downloadedLdaConf=s.downloadedLdaConf;
  downloadedDccConf=s.downloadedDccConf;
  downloadedDifConf=s.downloadedDifConf;
  downloadedAsicConf=s.downloadedAsicConf;
  locked=s.locked;
  comment=s.comment;
  nodb = s.nodb;
}

/**
Destructor
**/
State::~State(){
  if(lda_version!=nullptr)
    delete lda_version;
  if(dcc_version!=nullptr)
    delete dcc_version;
  if(dif_version!=nullptr)
    delete dif_version;
  if(asic_version!=nullptr)
    delete asic_version;  
  if(downloadedLdaConf!=nullptr)
    delete downloadedLdaConf;  
  if(downloadedDccConf!=nullptr)
    delete downloadedDccConf;  
  if(downloadedDifConf!=nullptr)
    delete downloadedDifConf;  
  if(downloadedAsicConf!=nullptr)
    delete downloadedAsicConf;  
}

State& State::operator=(const State& s){
  name = s.name;
  if(s.lda_version!=nullptr)
    lda_version = s.lda_version->clone();
  else
    lda_version=nullptr;
  if(s.dcc_version!=nullptr)
    dcc_version = s.dcc_version->clone();
  else
    dcc_version=nullptr;
  if(s.dif_version!=nullptr)
    dif_version = s.dif_version->clone();
  else
    dif_version=nullptr;
  if(s.asic_version!=nullptr)
    asic_version = s.asic_version->clone();
  else
    asic_version=nullptr;
  downloadedLdaConf=s.downloadedLdaConf;
  downloadedDccConf=s.downloadedDccConf;
  downloadedDifConf=s.downloadedDifConf;
  downloadedAsicConf=s.downloadedAsicConf;
  locked=s.locked;
  comment=s.comment;
  nodb=s.nodb;
  return *this;
}

string State::getName(){
  return name;
}

bool State::isLocked(){
  return locked;
}

Version* State::getLdaVersion(){
  if(lda_version==nullptr)
    return nullptr;
  return lda_version->clone();
}

Version* State::getDccVersion(){
  if(dcc_version==nullptr)
    return nullptr;
  return dcc_version->clone();
}

Version* State::getDifVersion(){
  if(dif_version==nullptr)
    return nullptr;
  return dif_version->clone();
}

Version* State::getAsicVersion(){
  if(asic_version==nullptr)
    return nullptr;
  return asic_version->clone();
}

void State::setName(string n){
  name = n;
}

void State::setComment(string c){
  comment = c;
}

State* State::download(string n) throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();  
  State *s = new State();
  try{
    Statement *stmt = db->createStatement("select name, ldaMajorId, ldaMinorId, LdaCreation, LdaVersionName, LdaLocked, dccMajorId, dccMinorId, dccCreation, dccVersionName, dccLocked, difMajorId, difMinorId, difCreation, difVersionName, difLocked, asicMajorId, asicMinorId, asicCreation, asicVersionName, asicLocked, state_comment, locked from STATE_DETAIL WHERE name=:n");
    stmt->setString(1,n);
    ResultSet *rs = stmt->executeQuery();
    if(rs->next()){
      s->name = rs->getString(1);
      s->lda_version = new Version(rs->getInt(2), rs->getInt(3), rs->getString(4), rs->getString(5), rs->getInt(6));
      s->dcc_version = new Version(rs->getInt(7), rs->getInt(8), rs->getString(9), rs->getString(10), rs->getInt(11));
      s->dif_version = new Version(rs->getInt(12), rs->getInt(13), rs->getString(14), rs->getString(15), rs->getInt(16));
      s->asic_version = new Version(rs->getInt(17), rs->getInt(18), rs->getString(19), rs->getString(20), rs->getInt(21));
      s->comment = rs->getString(22);
      s->locked = (rs->getInt(23)==1);
     }
    else{ 
      stmt->closeResultSet(rs);
      db->terminateStatement(stmt);
      ostringstream oss;
      oss<<"State "<<n<<" not found!";
      throw ILCException::Exception(oss.str());
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
  catch(ILCException::Exception e){
    throw e;
  }
  return s;
}

vector<int> State::getDifIDs(){
  if(nodb){
    throw ILCException::Exception("State::getDifIDs : This state has been created with no access to the database. It can not interact with the DB.");
  }
  DBConnection *db = DBConnection::getConnection();  
  vector<int> list;
  try{
    if(dif_version==nullptr)
      return list;
    Statement *stmt = db->createStatement("select DIF_ID FROM DIF_CONFIG where majorVersionId=:Maj and minorVersionId=0 ORDER BY DIF_ID");
    stmt->setInt(1,dif_version->getMajorId());
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      list.push_back(rs->getInt(1));
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
  catch(ILCException::Exception e){
    throw e;
  }
  return list;
}

vector<int> State::getAsicIDs(){
  if(nodb){
    throw ILCException::Exception("State::getAsicIDs : This state has been created with no access to the database. It can not interact with the DB.");
  }
  DBConnection *db = DBConnection::getConnection();  
  vector<int> list;
  try{
    if(asic_version==nullptr)
      return list;
    Statement *stmt = db->createStatement("select ASIC_ID FROM ASIC_CONFIG where majorVersionId=:Maj and minorVersionId=0 ORDER BY ASIC_ID");
    stmt->setInt(1,asic_version->getMajorId());
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      list.push_back(rs->getInt(1));
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
  catch(ILCException::Exception e){
    throw e;
  }
  return list;
}

vector<int> State::getDccIDs(){
  if(nodb){
    throw ILCException::Exception("State::getDccIDs : This state has been created with no access to the database. It can not interact with the DB.");
  }
  DBConnection *db = DBConnection::getConnection();  
  vector<int> list;
  try{
    if(dcc_version==nullptr)
      return list;
    Statement *stmt = db->createStatement("select DCC_ID FROM DCC_CONFIG where majorVersionId=:Maj and minorVersionId=0 ORDER BY DCC_ID");
    stmt->setInt(1,dcc_version->getMajorId());
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      list.push_back(rs->getInt(1));
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
  catch(ILCException::Exception e){
    throw e;
  }
  return list;
}

vector<string> State::getLdaIDs(){
  if(nodb){
    throw ILCException::Exception("State::getLdaIDs : This state has been created with no access to the database. It can not interact with the DB.");
  }
  DBConnection *db = DBConnection::getConnection();  
  vector<string> list;
  try{
    if(lda_version==nullptr)
      return list;
    Statement *stmt = db->createStatement("select LDA_ID FROM LDA_CONFIG where majorVersionId=:Maj and minorVersionId=0 ORDER BY LDA_ID");
    stmt->setInt(1,lda_version->getMajorId());
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      list.push_back(rs->getString(1));
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
  catch(ILCException::Exception e){
    throw e;
  }
  return list;
}

vector<int> State::getUsedDifIDs(){
  if(nodb){
    throw ILCException::Exception("State::getUsedDifIDs : This state has been created with no access to the database. It can not interact with the DB.");
  }
  DBConnection *db = DBConnection::getConnection();  
  vector<int> list;
  try{
    Statement *stmt = db->createStatement("select DISTINCT DIF_ID FROM ASIC_CONFIG, ASIC where ASIC.ID=ASIC_CONFIG.ASIC_ID AND ASIC_CONFIG.majorVersionId=:Maj and ASIC_CONFIG.minorVersionId=0 ORDER BY DIF_ID");
    stmt->setInt(1,asic_version->getMajorId());
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      list.push_back(rs->getInt(1));
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
  catch(ILCException::Exception e){
    throw e;
  }
  return list;
}

vector<int> State::getUsedDccIDs(){
  if(nodb){
    throw ILCException::Exception("State::getUsedDccIDs : This state has been created with no access to the database. It can not interact with the DB.");
  }
  DBConnection *db = DBConnection::getConnection();  
  vector<int> list;
  try{
    Statement *stmt = db->createStatement("select DISTINCT DCC_ID FROM DIF_CONFIG where majorVersionId=:Maj and minorVersionId=0 AND DCC_ID IS NOT nullptr ORDER BY DCC_ID");
    stmt->setInt(1,dif_version->getMajorId());
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      list.push_back(rs->getInt(1));
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
  catch(ILCException::Exception e){
    throw e;
  }
  return list;
}
vector<string> State::getUsedLdaIDs(){
  if(nodb){
    throw ILCException::Exception("State::getUsedLdaIDs : This state has been created with no access to the database. It can not interact with the DB.");
  }
  DBConnection *db = DBConnection::getConnection();  
  vector<string> list;
  try{
    Statement *stmt = db->createStatement("select DISTINCT LDA_ADDRESS FROM DIF_CONFIG where majorVersionId=:Maj and minorVersionId=0 AND LDA_ADDRESS IS NOT nullptr ORDER BY LDA_ADDRESS");
    stmt->setInt(1,dif_version->getMajorId());
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      list.push_back(rs->getString(1));
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
  catch(ILCException::Exception e){
    throw e;
  }
  return list;
}

vector<string> State::getUsedLdaByDccIDs(){
  if(nodb){
    throw ILCException::Exception("State::getUsedLdaByDccIDs : This state has been created with no access to the database. It can not interact with the DB.");
  }
  DBConnection *db = DBConnection::getConnection();  
  vector<string> list;
  try{
    Statement *stmt = db->createStatement("select DISTINCT LDA_ADDRESS FROM DCC, DCC_CONFIG where DCC.ID=DCC_CONFIG.DCC_ID AND majorVersionId=:Maj and minorVersionId=0 AND LDA_ADDRESS IS NOT nullptr ORDER BY LDA_ADDRESS");
    stmt->setInt(1,dcc_version->getMajorId());
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      list.push_back(rs->getString(1));
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
  catch(ILCException::Exception e){
    throw e;
  }

  return list;
}

vector<string> State::getUsedLdaChannelByDcc(){
  if(nodb){
    throw ILCException::Exception("State::getUsedLdaChannelByDcc : This state has been created with no access to the database. It can not interact with the DB.");
  }
  DBConnection *db = DBConnection::getConnection();  
  vector<string> list;
  try{
    Statement *stmt = db->createStatement("select D.ID, D.LDA_ADDRESS||'/'||D.LDA_CHANNEL, DC.majorversionid, DC.minorversionid FROM DCC D, DCC_CONFIG DC WHERE D.ID=DC.DCC_ID AND DC.majorversionid=:maj AND DC.minorversionid<=:min ORDER BY D.ID, DC.majorversionid, DC.minorversionid DESC");
    stmt->setInt(1,dcc_version->getMajorId());
    stmt->setInt(2,dcc_version->getMinorId());
    ResultSet *rs = stmt->executeQuery();
    string previousID = "";
    while(rs->next()){
      string id = rs->getString(1);
      if(previousID!=id){
	string n = rs->getString(2);
	list.push_back(n);
      }
      previousID=id;
    }
  }
  catch(ILCException::Exception e){
    throw e;
  }
  return list;
}

vector<string> State::getUsedLdaChannelByDif(){
  if(nodb){
    throw ILCException::Exception("State::getUsedLdaChannelByDif : This state has been created with no access to the database. It can not interact with the DB.");
  }
  DBConnection *db = DBConnection::getConnection();  
  vector<string> list;
  try{
    Statement *stmt = db->createStatement("select DIF_ID, LDA_ADDRESS||'/'||LDA_CHANNEL, majorversionid, minorversionid FROM DIF_CONFIG WHERE LDA_ADDRESS IS NOT nullptr AND majorversionid=:maj AND minorversionid<=:min ORDER BY DIF_ID, majorversionid, minorversionid DESC");
    stmt->setInt(1,dif_version->getMajorId());
    stmt->setInt(2,dif_version->getMinorId());
    ResultSet *rs = stmt->executeQuery();
    string previousID = "";
    while(rs->next()){
      string id = rs->getString(1);
      if(previousID!=id){
	string n = rs->getString(2);
	list.push_back(n);
      }
      previousID=id;
    }
  }
  catch(ILCException::Exception e){
    throw e;
  }
  return list;
}

LdaConfiguration* State::getLdaConfiguration() throw (ILCException::Exception){
  try{
    if(downloadedLdaConf==nullptr){
      downloadedLdaConf = new LdaConfiguration();
      downloadedLdaConf->download(lda_version->getMajorId(), lda_version->getMinorId());
    }
    return downloadedLdaConf;
  }
  catch(ILCException::Exception e){
    throw e;
  }
}

DccConfiguration* State::getDccConfiguration() throw (ILCException::Exception){
  try{
    if(downloadedDccConf==nullptr){
      downloadedDccConf = new DccConfiguration();
      downloadedDccConf->download(dcc_version->getMajorId(), dcc_version->getMinorId());
    }
    return downloadedDccConf;
  }
  
  catch(ILCException::Exception e){
    throw e;
  }
}

DifConfiguration* State::getDifConfiguration() throw (ILCException::Exception){
  try{
    if(downloadedDifConf==nullptr){
      downloadedDifConf = new DifConfiguration();
      downloadedDifConf->download(dif_version->getMajorId(), dif_version->getMinorId());
    }
    return downloadedDifConf;
  }
  catch(ILCException::Exception e){
    throw e;
  }
}

AsicConfiguration* State::getAsicConfiguration() throw (ILCException::Exception){
  try{
    if(downloadedAsicConf==nullptr){
      downloadedAsicConf = new AsicConfiguration();
      downloadedAsicConf->download(asic_version->getMajorId(), asic_version->getMinorId());
    }
    return downloadedAsicConf;
  }
  catch(ILCException::Exception e){
    throw e;
  }
}

bool State::isValid(){
  if(nodb){
    throw ILCException::Exception("State::isValid : This state has been created with no access to the database. It can not interact with the DB.");
  }
  //We want a valid version for each object, aven if they are none of them in the configuration
  if(lda_version==nullptr || dcc_version==nullptr || asic_version==nullptr || dif_version==nullptr){
    throw ILCException::Exception("You should have a valid version for LDA, DCC, DIF and ASIC! If there are no LDA or DCC, use the 0.0 version.");
  }
  //If we have asics without DIFs -> no need to go on
  if(asic_version!=nullptr && dif_version==nullptr){
    throw ILCException::Exception("There are no DIF configured in this state! You can not use ASICs in a state without the corresponding DIFs!");
  }

  vector<string> usedChannelsByDcc;
  vector<string> usedChannelsByDif;
  set<string> channels;
  pair<set<string>::iterator,bool> ret;

  if(dif_version!=nullptr){

    // Check that the DIFs needed by the configured ASICs are present
    vector<int> difList=getDifIDs();
    vector<int> neededDif=getUsedDifIDs();

    uint cpt=0;
    for(uint i=0;i<neededDif.size();i++){
      while(cpt<difList.size() && difList[cpt]!=neededDif[i]){
	cpt++;
      }
      if(cpt>=difList.size()){
	ostringstream oss;
	oss<<"The Dif "<<neededDif[i]<<" is needed by an ASIC but does not seem to be configured!";
	throw ILCException::Exception(oss.str());
      }
    }

    //Check that the DCCs needed by the configured DIFs are present
    vector<int> neededDcc=getUsedDccIDs();
    vector<int> dccList=getDccIDs();

    cpt=0;
    for(uint i=0;i<neededDcc.size();i++){
      while(cpt<dccList.size() && dccList[cpt]!=neededDcc[i]){
	cpt++;
      }
      if(cpt>=dccList.size()){
	ostringstream oss;
	oss<<"The Dcc "<<neededDcc[i]<<" is needed by a DIF but does not seem to be configured!";
	throw ILCException::Exception(oss.str());
      }
    }

    //Check that the LDAs needed by the configured DIFs are present
    vector<string> neededLda=getUsedLdaIDs();
    vector<string> ldaList=getLdaIDs();

    cpt=0;
    for(uint i=0;i<neededLda.size();i++){
      while(cpt<ldaList.size() && (ldaList[cpt].compare(neededLda[i])!=0)){
	cpt++;
      }
      if(cpt>=ldaList.size()){
	ostringstream oss;
	oss<<"The LDA "<<neededLda[i]<<" is needed by a DIF but does not seem to be configured!";
	throw ILCException::Exception(oss.str());
      }
    }
    
    /*
    // Check the LDA Channels used by the DIFs are unique
    usedChannelsByDif = getUsedLdaChannelByDif();
    for(uint i=0;i<usedChannelsByDif.size();i++){
      ret = channels.insert(usedChannelsByDif[i]);
      if (ret.second==false){
	ostringstream oss;
	oss<<"The LDA channel "<<*ret.first<<" is used by at least 2 DIF!";
	throw ILCException::Exception(oss.str());
      }
    }
    */
  }

  if(dcc_version!=nullptr){
    //Check that the LDAs needed by the configured DCCs are present
    vector<string> neededLda=getUsedLdaByDccIDs();
    vector<string> ldaList=getLdaIDs();
    
    uint cpt=0;
    for(uint i=0;i<neededLda.size();i++){
      while(cpt<ldaList.size() && (ldaList[cpt].compare(neededLda[i])!=0)){
	cpt++;
      }
      if(cpt>=ldaList.size()){
	ostringstream oss;
	oss<<"The LDA "<<neededLda[i]<<" is needed by a DCC but does not seem to be configured!";
	throw ILCException::Exception(oss.str());
      }
    }

    /*
    //Check that the DCC are using different LDA channels
    vector<string> usedChannels = getUsedLdaChannelByDcc();
    for(uint i=0;i<usedChannelsByDcc.size();i++){
      ret = channels.insert(usedChannelsByDcc[i]);
      if (ret.second==false){
	ostringstream oss;
	oss<<"The LDA channel "<<*ret.first<<" is used by at least 2 DCC or DIF!";
	throw ILCException::Exception(oss.str());
      }
    }
    */

  }
  return true;
}

void State::setAsicVersion(Version* v){
  if(nodb){
    throw ILCException::Exception("State::setAsicVersion : This state has been created with no access to the database. It can not interact with the DB.");
  }
  if(v!=nullptr)
    asic_version=v->clone();
  else
    asic_version=nullptr;

  // delete the cached version
  if(downloadedAsicConf!=nullptr){
    delete downloadedAsicConf; 
    downloadedAsicConf = nullptr;
  } 
}

void State::setDifVersion(Version* v){
  if(nodb){
    throw ILCException::Exception("State::setDifVersion : This state has been created with no access to the database. It can not interact with the DB.");
  }
  if(v!=nullptr)
    dif_version=v->clone();
  else
    dif_version=nullptr;

  // delete the cached version
  if(downloadedDifConf!=nullptr){
    delete downloadedDifConf; 
    downloadedDifConf = nullptr;
  }
}

void State::setDccVersion(Version* v){
  if(nodb){
    throw ILCException::Exception("State::setDccVersion : This state has been created with no access to the database. It can not interact with the DB.");
  }
  if(v!=nullptr)
    dcc_version=v->clone();
  else
    dcc_version=nullptr;

  // delete the cached version
  if(downloadedDccConf!=nullptr){
    delete downloadedDccConf; 
    downloadedDccConf = nullptr;
  }
}

void State::setLdaVersion(Version* v){
  if(nodb){
    throw ILCException::Exception("State::setLdaVersion : This state has been created with no access to the database. It can not interact with the DB.");
  }
  if(v!=nullptr)
    lda_version=v->clone();
  else
    lda_version=nullptr;
 
  // delete the cached version
  if(downloadedLdaConf!=nullptr){
    delete downloadedLdaConf; 
    downloadedLdaConf = nullptr;
  }
}

void State::lock(){
  if(nodb){
    throw ILCException::Exception("State::lock : This state has been created with no access to the database. It can not interact with the DB.");
  }

    DBConnection *db = DBConnection::getConnection();
    Statement *stmt = db->createStatement("UPDATE STATE SET locked=1 WHERE name=:n");
    stmt->setString(1, name);
    stmt->executeUpdate();
    db->terminateStatement(stmt);
    locked=true;

  
}

vector<State*> State::getStates(){
  vector<State*> v;
  DBConnection *db = DBConnection::getConnection();  
  try{
    Statement *stmt = db->createStatement("select name, ldaMajorId, ldaMinorId, LdaCreation, LdaVersionName, LdaLocked, dccMajorId, dccMinorId, dccCreation, dccVersionName, dccLocked, difMajorId, difMinorId, difCreation, difVersionName, difLocked, asicMajorId, asicMinorId, asicCreation, asicVersionName, asicLocked, state_comment, locked from STATE_DETAIL ORDER BY name");
    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      State *s = new State();
      s->name = rs->getString(1);
      s->lda_version = new Version(rs->getInt(2), rs->getInt(3), rs->getString(4), rs->getString(5), rs->getInt(6));
      s->dcc_version = new Version(rs->getInt(7), rs->getInt(8), rs->getString(9), rs->getString(10), rs->getInt(11));
      s->dif_version = new Version(rs->getInt(12), rs->getInt(13), rs->getString(14), rs->getString(15), rs->getInt(16));
      s->asic_version = new Version(rs->getInt(17), rs->getInt(18), rs->getString(19), rs->getString(20), rs->getInt(21));
      s->comment = rs->getString(22);
      s->locked = (rs->getInt(23)==1);
      v.push_back(s);
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
  }
  
  catch(ILCException::Exception e){
    throw e;
  }
  return v;

}

State* State::clone(){
  if(nodb){
    throw ILCException::Exception("State::clone : This state has been created with no access to the database. It can not interact with the DB.");
  }
  State* s = new State(name);
  s->lda_version=lda_version->clone();
  s->dcc_version=dcc_version->clone();
  s->dif_version=dif_version->clone();
  s->asic_version=asic_version->clone();
  s->comment=comment;
  s->locked=locked;
  return s;
}

bool State::toBeUpdated(){
  if(nodb){
    throw ILCException::Exception("State::toBeUpdated : This state has been created with no access to the database. It can not interact with the DB.");
  }
  if((downloadedLdaConf!=nullptr && downloadedLdaConf->isModified()) ||
     (downloadedDccConf!=nullptr && downloadedDccConf->isModified()) ||
     (downloadedDifConf!=nullptr && downloadedDifConf->isModified()) ||
     (downloadedAsicConf!=nullptr && downloadedAsicConf->isModified())){
    return true;
  }
  return false;
}

bool State::stateExist(string n){
  if(nodb){
    throw ILCException::Exception("State::stateExist : This state has been created with no access to the database. It can not interact with the DB.");
  }

    bool result = false;
    DBConnection *db = DBConnection::getConnection();  
    Statement *stmt = db->createStatement("select id from STATE WHERE name=:n");
    stmt->setString(1, n);
    ResultSet *rs = stmt->executeQuery();
    if(rs->next()){
      result = true;
    }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);
    return result;
 
  
}

vector<string> State::getUnusedStateNames() throw (ILCException::Exception){
  DBConnection *db = DBConnection::getConnection();  
  vector<string> names;
  try{
    Statement *stmt = db->createStatement("select name from State where id not in (select distinct state_id from setup_state) order by name");
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

void State::saveToXML(string fileName){
  if(nodb){
    throw ILCException::Exception("State::saveToXML : This state has been created with no access to the database. It can not interact with the DB.");
  }
  fstream f;
  bool exists=false;
  f.open(fileName.c_str(),ios::in);
  if( f.is_open() ){
    exists=true;
  }
  f.close();

  if(!exists){
    f.open(fileName.c_str(),ios::out);
    f<<"<CONFIG>"<<endl;
    f.close();

    this->getLdaConfiguration()->saveToXML(fileName);
    this->getDccConfiguration()->saveToXML(fileName);
    this->getDifConfiguration()->saveToXML(fileName);
    this->getAsicConfiguration()->saveToXML(fileName);

    f.open(fileName.c_str(),fstream::out | ios::app);
    f<<ConfigObject::cacheToXML();
    f<<"</CONFIG>"<<endl;
    f.close();
  }
}

/*State* State::createStateFromXML(string stateName, string fileName){
  size_t found=fileName.find_last_of("/\\");
  string name="";
  if(found==string::npos)
    name = fileName;//name of the file
  else
    name = fileName.substr(found+1);//name of the file

  State* s = new State(stateName);
  LdaConfiguration* lda_conf = new LdaConfiguration();
  lda_conf->addFromXML(fileName);
  if(lda_conf->size()>0){
    string n = name;
    while(Configuration::versionExist("LDA", n)){
      n = Version::incrementName(n);
    }
    Version* lda_version = lda_conf->upload(n);
    s->setLdaVersion(lda_version);
  }
  DccConfiguration* dcc_conf = new DccConfiguration();
  dcc_conf->addFromXML(fileName);
  if(dcc_conf->size()>0){
    string n = name;
    while(Configuration::versionExist("DCC", n)){
      n = Version::incrementName(n);
    }
    Version* dcc_version = dcc_conf->upload(n);
    s->setDccVersion(dcc_version);
  }
  DifConfiguration* dif_conf = new DifConfiguration();
  dif_conf->addFromXML(fileName);
  if(dif_conf->size()>0){
    string n = name;
    while(Configuration::versionExist("DIF", n)){
      n = Version::incrementName(n);
    }
    Version* dif_version = dif_conf->upload(n);
    s->setDifVersion(dif_version);
  }
  AsicConfiguration* asic_conf = new AsicConfiguration();
  asic_conf->addFromXML(fileName);
  if(asic_conf->size()>0){
    string n = name;
    while(Configuration::versionExist("ASIC", n)){
      n = Version::incrementName(n);
    }
    Version* asic_version = asic_conf->upload(n);
    s->setAsicVersion(asic_version);
  }

  return s;
}*/

/*State* State::createStateFromXML_NODB(string stateName, string fileName){
  //Load cache data from the XML file instead of DB
  ConfigObject::cacheFromXML(fileName);

  State* s = new State(stateName);
  s->nodb = true;
  s->downloadedLdaConf = new LdaConfiguration();
  s->downloadedLdaConf->addFromXML(fileName);
  s->downloadedDccConf = new DccConfiguration();
  s->downloadedDccConf->addFromXML(fileName);
  s->downloadedDifConf = new DifConfiguration();
  s->downloadedDifConf->addFromXML(fileName);
  s->downloadedAsicConf = new AsicConfiguration();
  s->downloadedAsicConf->addFromXML(fileName);
  
  //empty the cache
  ConfigObject::terminate();
  
  return s;
}*/

string State::checkLdaParams(string params) throw (ILCException::Exception){
  if(nodb){
    throw ILCException::Exception("State::checkLdaParams : This state has been created with no access to the database. It can not interact with the DB.");
  }
  string query = "select "+params+", count(*) from lda_config LDA, (select lda_id,majorversionid MAJ, max(minorversionid) MIN from lda_config where majorversionid=:maj AND minorversionid<=:min group by lda_id, majorversionid) V where LDA.lda_id=V.lda_id AND LDA.majorversionid=V.MAJ AND LDA.minorversionid=V.MIN group by ("+params+") order by "+params;
  return checkParams(query, lda_version, params);
}

string State::checkDccParams(string params) throw (ILCException::Exception){
  if(nodb){
    throw ILCException::Exception("State::checkDccParams : This state has been created with no access to the database. It can not interact with the DB.");
  }
  string query = "select "+params+", count(*) from dcc_config DCC, (select dcc_id,majorversionid MAJ, max(minorversionid) MIN from dcc_config where majorversionid=:maj AND minorversionid<=:min group by dcc_id, majorversionid) V, dcc D where D.id=DCC.dcc_id AND DCC.dcc_id=V.dcc_id AND DCC.majorversionid=V.MAJ AND DCC.minorversionid=V.MIN group by ("+params+") order by "+params;
  return checkParams(query, dcc_version, params);
}

string State::checkDifParams(string difType, string params) throw (ILCException::Exception){
  if(nodb){
    throw ILCException::Exception("State::checkDifParams : This state has been created with no access to the database. It can not interact with the DB.");
  }
  string query = "select "+params+", count(*) from "+difType+"_CONFIG ST, dif_config DIF, (select dif_id,majorversionid MAJ, max(minorversionid) MIN from dif_config where majorversionid=:maj AND minorversionid<=:min group by dif_id, majorversionid) V where ST."+difType+"_id=V.dif_id AND ST.majorversionid=V.MAJ AND ST.minorversionid=V.MIN AND DIF.dif_id=V.dif_id AND DIF.majorversionid=V.MAJ AND DIF.minorversionid=V.MIN group by ("+params+") order by "+params;
  return checkParams(query, dif_version, params);
}

string State::checkAsicParams(string asicType, string params) throw (ILCException::Exception){
  if(nodb){
    throw ILCException::Exception("State::checkAsicParams : This state has been created with no access to the database. It can not interact with the DB.");
  }
  string query = "select "+params+", count(*) from "+asicType+"_CONFIG ST, asic_config ASIC, (select asic_id,majorversionid MAJ, max(minorversionid) MIN from asic_config where majorversionid=:maj AND minorversionid<=:min group by asic_id, majorversionid) V, asic A where A.id=ASIC.asic_id AND ST."+asicType+"_id=V.asic_id AND ST.majorversionid=V.MAJ AND ST.minorversionid=V.MIN AND ASIC.asic_id=V.asic_id AND ASIC.majorversionid=V.MAJ AND ASIC.minorversionid=V.MIN group by ("+params+") order by "+params;
  return checkParams(query, asic_version, params);
}

string State::checkParams(string request, Version* v, string params) throw (ILCException::Exception){
  if(nodb){
    throw ILCException::Exception("State::checkParams : This state has been created with no access to the database. It can not interact with the DB.");
  }
  DBConnection *db = DBConnection::getConnection();  
  vector<string> parametersName;
  ostringstream result;
  std::transform(params.begin(), params.end(), params.begin(), ::toupper);
  if(params.find("INSERT")!=string::npos || 
     params.find("SELECT")!=string::npos || 
     params.find("UPDATE")!=string::npos || 
     params.find("DELETE")!=string::npos ||
     params.find(";")!=string::npos || 
     params.find("(")!=string::npos){
    ostringstream oss;
    oss<<"State::checkParams : bad format of params ("<<params<<")!";
    throw ILCException::Exception(oss.str());
  }

  size_t found;
  size_t previous=0;
  int count = 0;
  found=params.find_first_of(",");
  parametersName.push_back(params.substr(previous,found));
  while (found!=string::npos)
  {
    count++;
    previous=found;
    found=params.find_first_of(",",found+1);
    parametersName.push_back(params.substr(previous+1,found-previous-1));
  }
  count+=2;


    Statement *stmt = db->createStatement(request);
    stmt->setInt(1,v->getMajorId());
    stmt->setInt(2,v->getMinorId());

    for(uint j=0;j<parametersName.size();j++){
      result<<parametersName[j]<<"\t";
    }    
    result<<"NB"<<endl;

    ResultSet *rs = stmt->executeQuery();
    while(rs->next()){
      for(int i=1;i<count+1;i++){
	result<<rs->getString(i)<<"\t";
      }
      result<<endl;
     }
    stmt->closeResultSet(rs);
    db->terminateStatement(stmt);

 
  return result.str();
}
