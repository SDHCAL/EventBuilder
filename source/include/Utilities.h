#ifndef UTILITIES_H
#define UTILITIES_H
#include <iostream> 
#include"Colors.h"
#include "TH3.h"
#include "TMath.h"
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
void FillDelimiter(std::string &ToParse,int size,std::map<int ,std::vector<double> >&Delimiter)
{
    std::string delimiter_Dif = "|";
    std::string delimiter_Difs = "*";
    std::string delimiter_others=":";
    size_t pos = 0;
    std::string token;
    std::vector<std::string>a;
    bool findstar =ToParse.find(delimiter_Difs);
    if(findstar) {
        pos=ToParse.find(delimiter_Dif);
        a.push_back(ToParse.substr(0, pos));
        std::vector<double> tab(4);
        //int Dif=0;
        int j =0;
        size_t posi =0;
        while ((posi = a[0].find(delimiter_others)) != std::string::npos) {
            std::string token;
            token = a[0].substr(0, posi);
            tab[j]=atof(token.c_str());/*std::cout <<green<< token <<normal<< std::endl;*/



            a[0].erase(0, posi + delimiter_others.length());
            ++j;
        }
        tab[3]=atof(a[0].c_str());
        //std::cout<<green<<tab[3]<<normal<<std::endl;
        for(/*unsigned*/ int i=0; i<size; ++i) Delimiter[i+1]=tab;
    } else {
        while ((pos = ToParse.find(delimiter_Dif)) != std::string::npos) {
            token = ToParse.substr(0, pos);
            //std::cout << token << std::endl;
            a.push_back(token);
            ToParse.erase(0, pos + delimiter_Dif.length());
        }
        for(unsigned int i=0; i<a.size(); ++i) {
            std::vector<double> tab(4);
            int Dif=0;
            int j =0;
            size_t posi =0;
            while ((posi = a[i].find(delimiter_others)) != std::string::npos) {
                token = a[i].substr(0, posi);

                if (j==0) {
                    Dif=atof(token.c_str());/*std::cout <<blue<< token <<normal<< std::endl;*/
                } else {
                    tab[j-1]=atof(token.c_str());/*std::cout <<green<< token <<normal<< std::endl;*/
                }
                a[i].erase(0, posi + delimiter_others.length());
                ++j;
            }
            tab[3]=atof(a[i].c_str());
            //std::cout<<green<<tab[3]<<normal<<std::endl;
            Delimiter[Dif]=tab;
        }
    }
    if(a.size()==0)std::cout<<red<<"Warning:No Delimiters given "<<normal<<std::endl;
    if(a.size()!=(unsigned int)size&&a.size()!=0&&!findstar) {
        std::cout<<red<<"Error:Delimiters no well set ! "<<normal<<std::endl;
        std::exit(2);
    }

    //std::cout<<red<<a.size()<<normal<<std::endl;

    std::vector<std::string>word {"Imin : "," Imax : "," Jmin : "," Jmax : "};
    std::cout<<green<<"Delimiters"<<normal<<std::endl ;
    for(std::map<int,std::vector<double>>::iterator it=Delimiter.begin(); it!=Delimiter.end(); ++it) {

        std::cout<<green<<"Plane "<<it->first<<" : "<<normal;
        for(unsigned int i=0; i<it->second.size(); ++i) {
            std::cout<<green<<word[i]<<(it->second)[i]<<normal;
        }
        std::cout<<normal<<std::endl;
    }

}





Double_t transfer_function(const Double_t *x, const Double_t * /*param*/)
{
   if (*x <=0)return 0.00;
   else if(*x<=0.1) return 0.1;
   else if(*x<=0.2) return 0.2;
   else if(*x<=0.3) return 0.3;
   else if(*x<=0.4) return 0.4;
   else if(*x<=0.5) return 0.5;
   else if(*x<=0.6) return 0.6;
   else if(*x<=0.7) return 0.7;
   else if(*x<=0.8) return 0.8;
   else if(*x<=0.9) return 0.9;
   else return 1.0;
}

void Make_good_TH3(TH3* t)
{
 for(int i = 1; i <= t->GetNbinsX(); ++i)
    {
        	for(int j = 1; j <= t->GetNbinsY(); ++j)
        	{
            		for(int k = 1; k <= t->GetNbinsZ(); ++k)
            		{
                                
                		if(t->GetBinContent(i, j, k)!=0) t->SetBinContent(i, j, k,TMath::Log(t->GetBinContent(i, j, k)));
                                //std::cout<<hist->GetBinContent(i, j, k)<<std::endl;
            		}
        	}
    }
    long int val=0;
    for(int i = 1; i <= t->GetNbinsX(); ++i)
    {
        	for(int j = 1; j <= t->GetNbinsY(); ++j)
        	{
            		for(int k = 1; k <= t->GetNbinsZ(); ++k)
            		{
                                
                		val += t->GetBinContent(i, j, k);
            		}
        	}
    }
    for(int i = 1; i <= t->GetNbinsX(); ++i)
    {
        	for(int j = 1; j <= t->GetNbinsY(); ++j)
        	{
            		for(int k = 1; k <= t->GetNbinsZ(); ++k)
            		{
                                
                		if(t->GetBinContent(i, j, k)!=0) t->SetBinContent(i, j, k,t->GetBinContent(i, j, k)*1.0/val);
                                //std::cout<<hist->GetBinContent(i, j, k)<<std::endl;
            		}
        	}
    }
    double min=99999;
    double max=-1;
    for(int i = 1; i <= t->GetNbinsX(); ++i)
    {
        	for(int j = 1; j <= t->GetNbinsY(); ++j)
        	{
            		for(int k = 1; k <= t->GetNbinsZ(); ++k)
            		{
                                
                		if(t->GetBinContent(i, j, k)!=0) 
				{
					if(t->GetBinContent(i, j, k)<min)min=t->GetBinContent(i, j, k);
					if(t->GetBinContent(i, j, k)>max)max=t->GetBinContent(i, j, k);
				}
                                //std::cout<<hist->GetBinContent(i, j, k)<<std::endl;
            		}
        	}
    }
    for(int i = 1; i <= t->GetNbinsX(); ++i)
    {
        	for(int j = 1; j <= t->GetNbinsY(); ++j)
        	{
            		for(int k = 1; k <= t->GetNbinsZ(); ++k)
            		{
                                
                		if(t->GetBinContent(i, j, k)!=0) t->SetBinContent(i, j, k,t->GetBinContent(i, j, k)*1.0/max);
                                //std::cout<<hist->GetBinContent(i, j, k)<<std::endl;
            		}
        	}
    }
}

#endif
