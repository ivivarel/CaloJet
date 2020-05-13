#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TString.h>

#include <vector>

bool PrintHistogram(std::vector<TTree *> mytree /* array of histograms */, TString var, TString selection,  unsigned int nbins, float lowx, float highx);

void CompareSimulation(TString inputFileName1,TString inputFileName2)
{
  TFile * inputfile[2] = {0};
  inputfile[0] = TFile::Open(inputFileName1);
  inputfile[1] = TFile::Open(inputFileName2);

  std::vector<TTree *> mytree;
  mytree.push_back((TTree *) inputfile[0]->Get("MyTree"));
  mytree.push_back((TTree *) inputfile[1]->Get("MyTree"));

  TString selection = "j1t_theta > 0.7 && j1t_theta < 2.4 && j2t_theta > 0.7 && j2t_theta < 2.4 && j1t_E > 10 && j1s_E > 10 && j1c_E > 10 && j2t_E > 10 && j2s_E > 10 && j2c_E > 10";

  TString var = "higgsscint_m";

  if (!PrintHistogram(mytree,var,selection,50,100,150)){
    std::cout << "Cannot work with " << var << std::endl;
  }
  
  var = "higgscher_m";
  
  if (!PrintHistogram(mytree,var,selection,50,100,150)){
    std::cout << "Cannot work with " << var << std::endl;
  }
  
  var = "j1s_E";
  
  if (!PrintHistogram(mytree,var,selection,50,0,200)){
    std::cout << "Cannot work with " << var << std::endl;
  }
  
  var = "j1c_E";
  
  if (!PrintHistogram(mytree,var,selection,50,0,200)){
    std::cout << "Cannot work with " << var << std::endl;
  }

  var = "j1t_E";
  
  if (!PrintHistogram(mytree,var,selection,50,0,200)){
    std::cout << "Cannot work with " << var << std::endl;
  }

  var = "j2s_E";
  
  if (!PrintHistogram(mytree,var,selection,50,0,200)){
    std::cout << "Cannot work with " << var << std::endl;
  }

  var = "j2c_E";
  
  if (!PrintHistogram(mytree,var,selection,50,0,200)){
    std::cout << "Cannot work with " << var << std::endl;
  }

  var = "j2t_E";
  
  if (!PrintHistogram(mytree,var,selection,50,0,200)){
    std::cout << "Cannot work with " << var << std::endl;
  }

  var = "j1s_eta";
  
  if (!PrintHistogram(mytree,var,selection,60,-3,3)){
    std::cout << "Cannot work with " << var << std::endl;
  }

  var = "j1c_eta";
  
  if (!PrintHistogram(mytree,var,selection,50,-3,3)){
    std::cout << "Cannot work with " << var << std::endl;
  }

  var = "j1t_eta";
  
  if (!PrintHistogram(mytree,var,selection,50,-3,3)){
    std::cout << "Cannot work with " << var << std::endl;
  }




}

bool PrintHistogram(std::vector<TTree *> mytree /* array of histograms */, TString var, TString selection, unsigned int nbins, float lowx, float highx)
{

  if (mytree.size() != 2) {
    std::cout << "mytree.size() = " << mytree.size() << std::endl;
    return false;
  }
  
  TCanvas c1; 
  TH1F * h[2];
  TString hname;

  for (unsigned int i = 0; i < 2; ++i){
    hname = "h_" + var + "_";
    hname += i;
    h[i] = new TH1F(hname,"",nbins,lowx,highx);
    mytree[i]->Draw(var + " >> " + hname,selection);
    h[i]->SetLineColor(i+1);
    h[i]->SetXTitle(var);
    h[i]->Scale(1./h[i]->Integral(0,100000000));
  }

  h[0]->Draw("E");
  h[1]->Draw("HISTsame");
  c1.Print("SimulationComparison_IndependentGeneration/" + var + ".pdf");

  return true;
}
