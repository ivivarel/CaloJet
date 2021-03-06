//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 21 18:40:04 2013 by ROOT version 5.34/04
// from TTree MyTree/Ntuple
// found on file: StopSignal.root
//////////////////////////////////////////////////////////

#ifndef MyTree_h
#define MyTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class MyTree {
public :
   TTree           *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // ---- Declaration of leaf types ----
   Int_t nmuon ;
   Int_t nneu ;
   Double_t mjjr;
   Double_t mjjt;
   Double_t higgsrec_m;
   Double_t higgsrec_theta;
   Double_t higgsrec_E;
   Double_t higgsrec_phi;
   Double_t higgsrecem_m;
   Double_t higgsrecem_theta;
   Double_t higgsrecem_E;
   Double_t higgsrecem_phi;
   Double_t higgsHybErec_m;
   Double_t higgsHybErec_theta;
   Double_t higgsHybErec_E;
   Double_t higgsHybErec_phi;
   Double_t higgsHybDirRec_m;
   Double_t higgsHybDirRec_theta;
   Double_t higgsHybDirRec_E;
   Double_t higgsHybDirRec_phi;
   Double_t higgstruth_m;
   Double_t higgstruth_theta;
   Double_t higgstruth_E;
   Double_t higgstruth_phi;
   Double_t higgsscint_m;
   Double_t higgsscint_theta;
   Double_t higgsscint_E;
   Double_t higgsscint_phi;
   Double_t higgscher_m;
   Double_t higgscher_theta;
   Double_t higgscher_E;
   Double_t higgscher_phi;
   Double_t edep;
   Double_t muene_sci;
   Double_t muene_che;
   Double_t emcomp1;
   Double_t emcomp2;
   Double_t etotjr1;
   Double_t etotjr2;
   Double_t eleak;
   Double_t eleakn;
   Double_t drmmu;
   Double_t enumu;
   Double_t mnumu;
   Double_t emu;
   Double_t mbos_noc;

   Double_t j1t_E;
   Double_t j1t_pt;
   Double_t j1t_eta;
   Double_t j1t_phi;
   Double_t j1t_m;
   Double_t j1t_theta;
   Double_t j2t_E;
   Double_t j2t_pt;
   Double_t j2t_eta;
   Double_t j2t_phi;
   Double_t j2t_m;
   Double_t j2t_theta;
   Double_t j1r_E;
   Double_t j1r_pt;
   Double_t j1r_eta;
   Double_t j1r_phi;
   Double_t j1r_m;
   Double_t j1r_theta;
   Double_t j2r_E;
   Double_t j2r_pt;
   Double_t j2r_eta;
   Double_t j2r_phi;
   Double_t j2r_m;
   Double_t j2r_theta;
   Double_t j1s_E;
   Double_t j1s_pt;
   Double_t j1s_eta;
   Double_t j1s_phi;
   Double_t j1s_m;
   Double_t j1s_theta;
   Double_t j2s_E;
   Double_t j2s_pt;
   Double_t j2s_eta;
   Double_t j2s_phi;
   Double_t j2s_m;
   Double_t j2s_theta;
   Double_t j1c_E;
   Double_t j1c_pt;
   Double_t j1c_eta;
   Double_t j1c_phi;
   Double_t j1c_m;
   Double_t j1c_theta;
   Double_t j2c_E;
   Double_t j2c_pt;
   Double_t j2c_eta;
   Double_t j2c_phi;
   Double_t j2c_m;
   Double_t j2c_theta;
   Double_t j1s_eshare;
   Double_t j2s_eshare;
   Double_t j1c_eshare;
   Double_t j2c_eshare;
   Double_t closestT_DR_gam1;
   Double_t closestT_DR_gam2;
//

   TBranch *b_nmuon ;
   TBranch *b_nneu ;
   TBranch *b_higgsreco_m;
   TBranch *b_higgsreco_theta;
   TBranch *b_higgsreco_E;
  TBranch *b_higgsreco_phi;
  TBranch *b_higgsrecem_m;                                                
  TBranch *b_higgsrecem_theta;
  TBranch *b_higgsrecem_E;                                                                             
  TBranch *b_higgsrecem_phi;                                                                                      
  TBranch *b_higgsHybErec_m;                                                                                 
  TBranch *b_higgsHybErec_theta;                                                                             
  TBranch *b_higgsHybErec_E;                                                                              
  TBranch *b_higgsHybErec_phi;                                                                                  
  TBranch *b_higgsHybDirRec_m;                                                                                  
   TBranch *b_higgsHybDirRec_theta;                                                                               
   TBranch *b_higgsHybDirRec_E;                                                                                   
   TBranch *b_higgsHybDirRec_phi;  
   TBranch *b_higgstruth_m;
   TBranch *b_higgstruth_theta;
   TBranch *b_higgstruth_E;
   TBranch *b_higgstruth_phi;
   TBranch *b_higgsscint_m;
   TBranch *b_higgsscint_theta;
   TBranch *b_higgsscint_E;
   TBranch *b_higgsscint_phi;
   TBranch *b_higgscher_m;
   TBranch *b_higgscher_theta;
   TBranch *b_higgscher_E;
   TBranch *b_higgscher_phi;
   TBranch *b_mjjr;
   TBranch *b_mjjt;
   TBranch *b_edep;
   TBranch *b_muene_sci;
   TBranch *b_muene_che;
   TBranch *b_emcomp1;
   TBranch *b_emcomp2;
   TBranch *b_etotjr1;
   TBranch *b_etotjr2;
   TBranch *b_eleak;
   TBranch *b_eleakn;
   TBranch *b_drmmu;
   TBranch *b_enumu;
   TBranch *b_mnumu;
   TBranch *b_emu;
   TBranch *b_mbos_noc;

   TBranch *b_j1t_E;
   TBranch *b_j1t_pt;
   TBranch *b_j1t_eta;
   TBranch *b_j1t_phi;
   TBranch *b_j1t_m;
   TBranch *b_j1t_theta;
   TBranch *b_j2t_E;
   TBranch *b_j2t_pt;
   TBranch *b_j2t_eta;
   TBranch *b_j2t_phi;
   TBranch *b_j2t_m;
   TBranch *b_j2t_theta;
   TBranch *b_j1r_E;
   TBranch *b_j1r_pt;
   TBranch *b_j1r_eta;
   TBranch *b_j1r_phi;
   TBranch *b_j1r_m;
   TBranch *b_j1r_theta;
   TBranch *b_j2r_E;
   TBranch *b_j2r_pt;
   TBranch *b_j2r_eta;
   TBranch *b_j2r_phi;
   TBranch *b_j2r_m;
   TBranch *b_j2r_theta;
   TBranch *b_j1s_E;
   TBranch *b_j1s_pt;
   TBranch *b_j1s_eta;
   TBranch *b_j1s_phi;
   TBranch *b_j1s_m;
   TBranch *b_j1s_theta;
   TBranch *b_j2s_E;
   TBranch *b_j2s_pt;
   TBranch *b_j2s_eta;
   TBranch *b_j2s_phi;
   TBranch *b_j2s_m;
   TBranch *b_j2s_theta;
   TBranch *b_j1c_E;
   TBranch *b_j1c_pt;
   TBranch *b_j1c_eta;
   TBranch *b_j1c_phi;
   TBranch *b_j1c_m;
   TBranch *b_j1c_theta;
   TBranch *b_j2c_E;
   TBranch *b_j2c_pt;
   TBranch *b_j2c_eta;
   TBranch *b_j2c_phi;
   TBranch *b_j2c_m;
   TBranch *b_j2c_theta;
   TBranch *b_closestT_DR_gam1;
   TBranch *b_closestT_DR_gam2;

   MyTree();
   virtual ~MyTree();
   void     Init(); 
   virtual void     Write();
   void     Fill();
   void     Reset();
   
};

#endif

#ifdef MyTree_cxx
MyTree::MyTree() 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   //if (tree == 0) {
   /*
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("StopSignal.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("StopSignal.root");
      }
      f->GetObject("MyTree",tree);
   */
   //}
   Init();
}

MyTree::~MyTree()
{
}

#endif // #ifdef MyTree_cxx
