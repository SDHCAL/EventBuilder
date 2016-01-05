#ifndef _SETUP_H_
#define _SETUP_H_

#include <iostream>
#include <sstream>
#include <string>
#include <set>

#include "../db/DBConnection.h"
#include "State.h"

using namespace std;

/**
  A Setup is a set of States. While each state contains the description of a subdetector configuration, the setup sumup the whole detector configuration.
**/
class Setup{

 private :
  string name;
  bool locked;
  vector<State*> states;

  bool isValid();
  bool setupExist();

 public:
  /**
     Constructor
  **/
  Setup(string n);
  ~Setup();
  Setup(Setup& s);

  /**
     Add a State to the setup (the State is cloned)
  **/
  void addState(State* s);
    /**
     Get the list of states associated to this setup
     @return The vector of states. This is NOT a clone! Any modification on this vector will modify the setup. You should not delete it, this vector will be destroy with the setup object.
  **/
  vector<State*> getStates();

  /**
     Get the name of the setup
     @return The name of the setup
  **/
  string getName();

  /**
   Lock the setup so that it can not be modified
  **/
  void lock();
  /**
     Check if the setup is locked
     @return True if the setup is locked
  **/
  bool isLocked();
  
  /**
     Retrieve a setups from the database
     @param n The name of the setup
     @return A pointer on the setup object (you have to delete it)
  **/
  static Setup* getSetup(string n) throw (ILCException::Exception);

  /**
     Retrieve all the setup names from the database
     @return A vector of setup's name
  **/
  static vector<string> getSetupNames() throw (ILCException::Exception);

  /**
     Retrieve the names of the setups not used in a run
     @return A vector of setup's name
  **/
  static vector<string> getUnusedSetupNames() throw (ILCException::Exception);
  
  /**
     Upload the current setup to the database (update it if already exists)
  **/
  void upload() throw (ILCException::Exception);

  /**
     Get the setup used for a given run
     @param runN The run number
     @return A pointer on the setup used for the run
  **/
  static Setup* getSetupFromRun(int runN) throw (ILCException::Exception);

};
#endif
