#ifndef INTRO_H
#define INTRO_H
#include <iostream>
#include "Version.h"
#include "Colors.h"
#include <unistd.h>
#include <cstdlib>
void inline Intro()
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<red<<" ______               _   ____        _ _     _            "<<normal<<std::endl;          
  std::cout<<red<<"|  ____|             | | |  _ \\      (_) |   | |           "<<normal<<std::endl;  
  std::cout<<red<<"| |____   _____ _ __ | |_| |_) |_   _ _| | __| | ___ _ __  "<<normal<<std::endl;
  std::cout<<red<<"|  __\\ \\ / / _ \\ '_ \\| __|  _ <| | | | | |/ _` |/ _ \\ '__| "<<normal<<std::endl;
  std::cout<<red<<"| |___\\ V /  __/ | | | |_| |_) | |_| | | | (_| |  __/ |    "<<normal<<std::endl;
  std::cout<<red<<"|______\\_/ \\___|_| |_|\\__|____/ \\__,_|_|_|\\__,_|\\___|_| ";
  std::cout<<green<<"v "<<Test_VERSION_MAJOR<<"."<<Test_VERSION_MINOR<<normal<<std::endl;  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;
  sleep(1);                                        
}
#endif
