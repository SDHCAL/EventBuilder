#include <iostream>
#include <string>
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <TApplication.h>
int main(int argc, char* argv[])
{
  TApplication* rootapp = new TApplication("example",&argc, argv);
  TCanvas* canvas = new TCanvas("canvas");
  TMultiGraph *c2 = new TMultiGraph();
  //std::cout<<"Titre : ";
  //std::string title="";
  //std::cin>>title;
  //const char* bbb=title.c_str();
  c2->SetTitle(argv[1]);
  TLegend* leg = new TLegend(1.0,0.8,0.9,0.7);
  
  TGraphErrors *c1 = new TGraphErrors("To_Plot.txt", "%lg %lg %lg");
  c1->SetLineColor(1);
  leg->AddEntry(c1, "Plane 1", "l");
 

   
  c2->Add(c1);
  
  
	
  std::string pp="%*lg %*lg %*lg %*lg";
  std::string dd="%*lg %*lg %*lg %*lg";
  std::cout<<argv[2]<<"  "<<argv[3]<<std::endl;
  for(int i=0;i<atoi(argv[2])-1;i++)
  {
	//std::cout<<i<<std::endl;
	std::string v="%lg "+dd+" %lg %lg ";
	//std::cout<<v<<std::endl;
        const char* bb=v.c_str();
	TGraphErrors *c1 = new TGraphErrors("To_Plot.txt", bb);
	
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
  rootapp->Run();

  return 0;
}

