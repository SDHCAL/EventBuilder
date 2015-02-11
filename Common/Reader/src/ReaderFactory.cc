#include "../include/ReaderFactory.h"
#include"../include/XMLReader.h"
#include"../include/TXTReader.h"

ReaderFactory::ReaderFactory()
{
	RegisterReader("XMLReader",new XMLReader());
  RegisterReader("TXTReader",new TXTReader());
}
