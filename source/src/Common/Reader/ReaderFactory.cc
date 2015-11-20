#include "Reader/ReaderFactory.h"
#include"Reader/XMLReader.h"
#include"Reader/TXTReader.h"
#include"Reader/XMLReaderConfig.h"

ReaderFactory::ReaderFactory()
{
	RegisterReader("XMLReader",new XMLReader());
  RegisterReader("TXTReader",new TXTReader());
  RegisterReader("XMLReaderConfig",new XMLReaderConfig());
}
