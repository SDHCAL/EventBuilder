// dump content of files
// done from tutorials copyFiles.C by Author: Rene Brun

#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

  //Example of script showing how to copy all objects (including directories)
  //from a source file.
  //For each input file, a new directory is created in the current directory 
  //with the name of the source file.
  //After execution of:
  // root > .x copyFiles.C
  //the file result.root will contain 4 subdirectories:
  // "tot100.root", "hsimple.root", "hs1.root","hs2.root"
      
void DumpDir(TDirectory *source) {
   //copy all objects and subdirs of directory source as a subdir of the current directory   
   source->ls();
   //loop on all entries of this directory
   TKey *key;
   TIter nextkey(source->GetListOfKeys());
   while ((key = (TKey*)nextkey())) {
      const char *classname = key->GetClassName();
      TClass *cl = gROOT->GetClass(classname);
      if (!cl) continue;
      if (cl->InheritsFrom(TDirectory::Class())) {
         source->cd(key->GetName());
         TDirectory *subdir = gDirectory;
         DumpDir(subdir);
      } else if (cl->InheritsFrom(TTree::Class())) {
         TTree *T = (TTree*)source->Get(key->GetName());
         T->Print();
	 T->SetScanField(0);
	 T->Scan("*");
      } else {
         source->cd();
         TObject *obj = key->ReadObj();
         obj->Print("ALL");
         delete obj;
     }
   }
}

void DumpFile(const char *fname) {
   TFile *f = TFile::Open(fname);
   if (!f || f->IsZombie()) {
      printf("Cannot dump file: %s\n",fname);
      return;
   }
   DumpDir(f);
   delete f;
}  
void dumpFiles() {
   //main function dumping 2 files
   DumpFile("Results_728770.root");
   DumpFile("Results_Analysis728770.root");   
}
