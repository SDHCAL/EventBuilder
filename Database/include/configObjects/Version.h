#ifndef _VERSION_H_
#define _VERSION_H_

#include<occi.h>
#include <string>
#include "db/DBConnection.h"

using namespace std;

/**
   Class used to manage versions
**/
class Version{

 private:
  int majorId;
  int minorId;
  string name;
  string creationDate;
  bool locked;

 public:

  /** 
      Constructor
      @param maj MajorVersionId
      @param min MinorVersionId
  **/
  Version(int maj, int min);

  /** 
      Constructor
      @param maj MajorVersionId
      @param min MinorVersionId
      @param n Name of the version
      @param l The version is locked/unlocked
  **/
  Version(int maj, int min, string n, bool l);
  /** 
      Constructor
      @param maj MajorVersionId
      @param min MinorVersionId
      @param d Creation time of the version
      @param n Name of the version
      @param l The version is locked/unlocked
  **/
  Version(int maj, int min, string d, string n, bool l);
  ~Version();
  /** Selector **/
  int getMajorId();
  /** Selector **/
  int getMinorId();
  /** Selector **/
  string getName();
  /** Selector **/
  string getDate();
  /** Selector **/
  bool isLocked();
  /** 
      Actualy create the version in the database (you should not have to use this method).
  **/
  void create(string table);


  /**
     Lock this version in the database so that no modification can be made (you should not have to use this method)
  **/
  void lock(string table);

  /**
     Delete all pointers contained in the vector and clear the vector
     @param v The vector to clear
  **/
  static void clearVersionVector(vector<Version*> v);

  /**
     Return a clone of the object
  **/
  Version* clone();

  /**
     Increment the number at the end of the name ('_X'), add one ('_1') if none is existing
     @param n The initial name
     @return The incremented name
  **/
  static string incrementName(string n);

  /**
     Check if this version is the most up-to-date one in this major line
     @param table The name of the table to check
     @return True if there is no following minor version
  **/
  bool isMostRecent(string table);

};
#endif
