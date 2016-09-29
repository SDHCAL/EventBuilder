#include "ReaderFactory.h"
#include"XMLReader.h"
#include"TXTReader.h"
#include"XMLReaderConfig.h"
#include "XMLReaderElog.h"

ReaderFactory::ReaderFactory()
{
  RegisterReader("XMLReader",new XMLReader());
  RegisterReader("TXTReader",new TXTReader());
  RegisterReader("XMLReaderConfig",new XMLReaderConfig());
  RegisterReader("XMLReaderElog",new XMLReaderElog());
}
