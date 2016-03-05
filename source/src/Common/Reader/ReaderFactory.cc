#include "Reader/ReaderFactory.h"
#include"Reader/XMLReader.h"
#include"Reader/TXTReader.h"
#include"Reader/XMLReaderConfig.h"
#include "Reader/XMLReaderElog.h"

ReaderFactory::ReaderFactory()
{
  RegisterReader("XMLReader",new XMLReader());
  RegisterReader("TXTReader",new TXTReader());
  RegisterReader("XMLReaderConfig",new XMLReaderConfig());
  RegisterReader("XMLReaderElog",new XMLReaderElog());
}
