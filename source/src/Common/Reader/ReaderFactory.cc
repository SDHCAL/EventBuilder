#include "Reader/ReaderFactory.h"
#include"Reader/XMLReader.h"
#include"Reader/TXTReader.h"

ReaderFactory::ReaderFactory()
{
	RegisterReader("XMLReader",new XMLReader());
  RegisterReader("TXTReader",new TXTReader());
}
