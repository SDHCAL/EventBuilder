#!/bin/bash

cd /home/grenier/Travail/ilcsoft/v01-17-06/v01-17-06/
source init_ilcsoft.sh

cd EventBuilder
export MARLIN_DLL=${PWD}/lib/libAnalysis.so

cd donneesFrancois

echo "running Marlin"

Marlin --global.LCIOInputFiles="DHCAL_Trivent_000000_I0.slcio" --global.MaxRecordNumber=0 xml/Analysis.xml > output.txt

echo "dumping ROOT files"

root -b -q dumpFiles.C > rootfiledump.txt

mv *.txt output_txt


diff output_txtref output_txt



 
