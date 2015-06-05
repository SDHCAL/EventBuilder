#include <iostream>
#include <string>
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <TApplication.h>
int main(int argc, char* argv[])
{
   //TApplication* rootapp = new TApplication("example",&argc, argv);
  TCanvas* canvas = new TCanvas("canvas");
  TMultiGraph *c2 = new TMultiGraph();
  TMultiGraph *c3 = new TMultiGraph();
  //std::cout<<"Titre : ";
  //std::string title="";
  //std::cin>>title;
  //const char* bbb=title.c_str();
  c2->SetTitle(argv[1]);
  TLegend* leg = new TLegend(0.8, 0.6, 0.9, 0.9);
  
  TGraphErrors *c1 = new TGraphErrors("Results.txt", "%lg %lg %lg");
  c1->SetLineColor(1);
  leg->AddEntry(c1, "Plane 1", "l");
 

   
  c2->Add(c1);
  
  
	
  std::string pp="%*lg %*lg %*lg %*lg";
  std::string dd="%*lg %*lg %*lg %*lg";
  for(int i=0;i<6;i++)
  {
	//std::cout<<i<<std::endl;
	std::string v="%lg "+dd+" %lg %lg ";
	//std::cout<<v<<std::endl;
        const char* bb=v.c_str();
	TGraphErrors *c1 = new TGraphErrors("Results.txt", bb);
	
        dd+=pp;
	c1->SetLineColor(i+2);
        c2->Add(c1);
        std::string lll = "Plane "+std::to_string(i+2);
        const char* llll=lll.c_str();
	leg->AddEntry(c1,llll , "l");
      
  }
  
  c2->Draw("APELSAME");
  leg->Draw();
   canvas->Print("plots.pdf");
//rootapp->Run();

  return 0;
}

