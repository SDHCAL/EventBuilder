#ifndef _DBINIT_H_
#define _DBINIT_H_

#include "DBConnection.h"
#include "configObjects/ConfigObject.h"

class DBInit{
 public:
  static void init();
  static void terminate();
 private:
  /** Brief retreive login, password and path from env. variable CONFDB
   * \param login - on the database
   * \param passwd - on the database
   * \param path - on database
   */
  static void getDbConfiguration(std::string &login, std::string &passwd, std::string &path);
};

#endif
