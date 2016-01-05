#ifndef _STATE_H_
#define _STATE_H_

#include <iostream>
#include <sstream>
#include <string>
#include <set>

#include "Version.h"
#include "../db/DBConnection.h"
#include "LdaConfiguration.h"
#include "DccConfiguration.h"
#include "DifConfiguration.h"
#include "AsicConfiguration.h"

using namespace std;

/**
   Class Managing a State
**/

class State{

 private:
  string name;
  Version* lda_version;
  Version* dcc_version;
  Version* dif_version;
  Version* asic_version;

  // When a user ask for a configuration, we keep a pointer on this new object for later use
  // This is to avoid multiple downloads and to keep a link on the configuration the user is possibly modifing
  // We will  delete these pointers when the state is deleted, so the user should not do it...
  LdaConfiguration* downloadedLdaConf;
  DccConfiguration* downloadedDccConf;
  DifConfiguration* downloadedDifConf;
  AsicConfiguration* downloadedAsicConf;

  string comment;
  bool locked;

  bool nodb;

  bool isValid();
  string checkParams(string request, Version* v, string params) throw (ILCException::Exception);

 public:
  /**
     Constructor
  **/
  State();
  /**
     Constructor
     @param n The name of the state
  **/
  State(string n);
  /**
     Copy Constructor
  **/
  State(const State& s);

  State& operator=(const State& s);
  ~State();

  /**
     Change the name of the state
     @param n The new name
  **/
  void setName(string n);
  /**
     Change the comment of the state
     @param n The new comment
  **/
  void setComment(string c);
 /**
     Change the LDA version of the state (the parameter is cloned!)
     Previously downloaded LdaConfigurations are set to NULL!
     @param v The new version
  **/
  void setLdaVersion(Version* v);
 /**
     Change the DCC version of the state (the parameter is cloned!)
     Previously downloaded DccConfigurations are set to NULL!
     @param v The new version
  **/
  void setDccVersion(Version* v);
 /**
     Change the DIF version of the state (the parameter is cloned!)
     Previously downloaded DifConfigurations are set to NULL!
     @param v The new version
  **/
  void setDifVersion(Version* v);
 /**
     Change the ASIC version of the state (the parameter is cloned!)
     Previously downloaded AsicConfigurations are set to NULL!
     @param v The new version
  **/
  void setAsicVersion(Version* v);

  string getName();
  
  /**
     Get the LDA version of the state
     @return A clone of the Version object (you must free the returned pointer once used!)
  **/
  Version* getLdaVersion();
  /**
     Get the DCC version of the state
     @return A clone of the Version object (you must free the returned pointer once used!)
  **/
  Version* getDccVersion();
  /**
     Get the DIF version of the state
     @return A clone of the Version object (you must free the returned pointer once used!)
  **/
  Version* getDifVersion();
  /**
     Get the ASIC version of the state
     @return A clone of the Version object (you must free the returned pointer once used!)
  **/
  Version* getAsicVersion();

  bool isLocked();

  /**
     Check the validity of the state and upload it to the database. May throw Exception.
  **/
  void upload() throw (ILCException::Exception);

  /**
    Download a state from the database. May throw Exception.
    @param n The name of the state to download
    @return A State object containing all the informations from the database (you have to delete it!).
  **/
  static State* download(string n) throw (ILCException::Exception);


  vector<int> getDifIDs();
  vector<int> getAsicIDs();
  vector<int> getDccIDs();
  vector<string> getLdaIDs();
  vector<int> getUsedDifIDs();
  vector<int> getUsedDccIDs();
  vector<string> getUsedLdaIDs();
  vector<string> getUsedLdaByDccIDs();
  vector<string> getUsedLdaChannelByDif();
  vector<string> getUsedLdaChannelByDcc();

  /**
     Download the LDA Configuration associated with the current state from the database. May throw Exception
     @return A LdaConfiguration pointer. This pointer will be automaticaly deleted in the State destructor.
  **/
  LdaConfiguration* getLdaConfiguration() throw (ILCException::Exception);
  /**
     Download the DCC Configuration associated with the current state from the database. May throw Exception
     @return A DccConfiguration pointer. This pointer will be automaticaly deleted in the State destructor.
  **/
  DccConfiguration* getDccConfiguration() throw (ILCException::Exception);
  /**
     Download the DIF Configuration associated with the current state from the database. May throw Exception
     @return A DifConfiguration pointer.  This pointer will be automaticaly deleted in the State destructor.
  **/
  DifConfiguration* getDifConfiguration() throw (ILCException::Exception);
  /**
     Download the ASIC Configuration associated with the current state from the database. May throw Exception
     @return A AsicConfiguration pointer. This pointer will be automaticaly deleted in the State destructor.
  **/
  AsicConfiguration* getAsicConfiguration() throw (ILCException::Exception);

  /**
     Lock the state so that no further modification can be made.
  **/
  void lock();

  /**
     Retrieve all the existing states in the database
     @return A vector containing all the existing states in the database. You have to delete all this pointers!
  **/
  vector<State*> getStates();

  /**
     Retrieve the names of the states not used in a setup
     @return A vector of setup's name
  **/
  static vector<string> getUnusedStateNames() throw (ILCException::Exception);

  /**
     Return a new State which is a copy of the current one
  **/
  State* clone();

  /**
     Check if the state does exist
     @param n The name of the state
     @return True if the state already exists in the DB
  **/
  bool stateExist(string n);

  /**
     Check if the state has been modified
  **/
  bool toBeUpdated();

  /**
     Export the configuration to a XML file
     @param fileName The name of the file to write the configuration (should not exist before the call).
  **/
  void saveToXML(string fileName);

  /**
     Get the values of the given parameters in the current state
     @param params A String containing the list of parameters separated with a coma (ex : 'FIRMWARE, HARDWARE')
     @return A string containing the values of the parameters + the number of LDA having such values
   **/
  string checkLdaParams(string params) throw (ILCException::Exception);
  /**
     Get the values of the given parameters in the current state
     @param params A String containing the list of parameters separated with a coma (ex : 'FIRMWARE, HARDWARE')
     @return A string containing the values of the parameters + the number of DCC having such values
   **/
  string checkDccParams(string params) throw (ILCException::Exception);
  /**
     Get the values of the given parameters in the current state
     @param difType The type of DIF
     @param params A String containing the list of parameters separated with a coma (ex : 'FIRMWARE, HARDWARE')
     @return A string containing the values of the parameters + the number of DIF having such values
   **/
  string checkDifParams(string difType, string params) throw (ILCException::Exception);
  /**
     Get the values of the given parameters in the current state
     @param asicType The type of ASIC
     @param params A String containing the list of parameters separated with a coma (ex : 'FIRMWARE, HARDWARE')
     @return A string containing the values of the parameters + the number of ASIC having such values
   **/
  string checkAsicParams(string asicType, string params) throw (ILCException::Exception);
};

#endif
