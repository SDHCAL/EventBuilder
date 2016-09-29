#!/bin/bash
CreateBuildDir()
{
  if [ -d build ]; then
	  rm -r build
  fi
  mkdir build
}

CreateMakefile()
{
    echo "============================================================"
    echo "=                 Build Makefile for compilation           ="
    echo "============================================================"    
    echo "  Create a new Makefile  .. "
    cmake ..
    ls -lthr  
}

Compile()
{
    echo "============================================================"
    echo "=                      Compilation                         ="
    echo "============================================================"
    make
    make install
}

#The script
CreateBuildDir
cd build
CreateMakefile
Compile
echo " Compilation done"
if [ -d build ]; then
	  rm -r build
fi
