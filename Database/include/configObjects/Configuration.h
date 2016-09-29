#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_

#include <vector>
#include <algorithm>
#include <fstream>
#include "ConfigObject.h"
#include "Version.h"

/**
   Abstract super class for configuration
**/
class Configuration{

 protected:
  string objectName;
  vector<ConfigObject*> list;
  bool modified;

  Version* createMajorVersion(string n) const throw (ILCException::Exception);
  Version* createMinorVersion(int m, string n) const throw (ILCException::Exception);

  /**
     Get the list of unused versions
     @param object_name The name of the configured object
     @return A vector containing all the unused versions. Once used, you should free the memory using the Version::clearVersionVector(vector<Version*> v) method.
  **/
  static vector<Version*> getUnusedVersions(string object_name) throw (ILCException::Exception);

 public:
  Configuration();
  virtual ~Configuration();
  
  /**
     Add an object to the configuration
     @param l A pointer to the object to add
  **/
  void add(ConfigObject* l);

  /**
     Get the vector actualy containing the configurations
     @return A vector containing all the configured objects (they are NOT clones!)
  **/
  vector<ConfigObject*>& getVector();

  /**
     Return the number of objects in the configuration
  **/
  int size();

  /**
     Delete all contained objects and clear the vector of objects
  **/
  void clear(); 
  /**
     Load the given version
     @param maj The majorVersionID
     @param min The MinorVersionID
  **/
  virtual void download(int maj, int min) = 0;
  
 /**
     Get all existing versions
     @param objectName Name of the object for which you want the existing versions (can be LDA, DCC, DIF, ...)
     @return A vector containing all the existing versions. Once used, you should free the memory using the Version::clearVersionVector(vector<Version*> v) method.
  **/
  static vector<Version*> getVersions(string objectName) throw (ILCException::Exception);

 /**
     Get an existing versions
     @param objectName Name of the object for which you want the existing version (LDA, DCC, DIF, ASIC, ...)
     @param majorId The version major ID
     @param minorId The version minor ID
     @return A pointer on a Version. Once used, you should free the memory with a delete on the pointer.
  **/
  static Version* getVersion(string objectName, int majorId, int minorId) throw (ILCException::Exception);


  /**
     Check if the version name is alreasy used
     @param objectName Name of the object for which you want the existing version (LDA, DCC, DIF, ASIC, ...)
     @param n The name of the version
     @return True if the name is already used by a version
  **/
  static bool versionExist(string objectName, string n) throw (ILCException::Exception);

  /**
     Flag the configuration as modified -> a new version configuration will be created in case of recursive upload of the containing state
     @param b True if the configuration was modified
  **/
  void setModified(bool b);
  
  /**
     Check the modified flag
     @return True if the configuration was flaged as modified
  **/
  bool isModified();

  /**
     Check if one of the objects was flagged as modified -> the user is preparing a minor update
  **/
  bool minorUpdateNeeded();

  /**
     Get a XML version of the configuration
     @return A vector<string> containing one XML description per object
  **/
  vector<string> toXML();

  /**
     Export the configuration to a XML file
     @param fileName The name of the file to write the configuration. If the file exists, the configuration is appended to the end. If the file does not exist, it is created.
  **/
  void saveToXML(string fileName);

};
#endif
