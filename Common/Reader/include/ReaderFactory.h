#ifndef READERFACTORY_H
#define READERFACTORY_H
#include"Reader.h"
#include<map>
#include<string>
class ReaderFactory
{
  public:
  ReaderFactory();
  ~ReaderFactory()
  {
    for(std::map<std::string, Reader *>::iterator iter = readerMap.begin(), endIter = readerMap.end() ; endIter != iter ; ++iter)delete iter->second;
    readerMap.clear();
  }
  Reader* CreateReader(const std::string &readerType) const
  {
     std::map<std::string, Reader*>::const_iterator findIter = readerMap.find(readerType);
     if(readerMap.end() == findIter)return NULL;
     else return findIter->second->Clone();
  }
  void RegisterReader(const std::string &readerType, Reader *pReader)
  {
    std::map<std::string, Reader*>::const_iterator findIter = readerMap.find(readerType);
    if(findIter != readerMap.end())return;
    readerMap[readerType] = pReader;
  }
  protected:
  std::map<std::string, Reader*> readerMap;
};
#endif
