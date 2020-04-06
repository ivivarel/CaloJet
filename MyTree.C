#define MyTree_cxx
#include "MyTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>

void MyTree::Write()
{
  fChain->Write();
}

void MyTree::Fill(){
  fChain->Fill();
 }


void MyTree::Init()
{
   fChain = new TTree("MyTree","Bonsai_Tree"); 
   
   fChain->Branch("nmuon ",&nmuon ,"nmuon/I");
   fChain->Branch("nneu ",&nneu ,"nneu/I");
   fChain->Branch("higgsrec_m",&higgsrec_m,"higgsrec_m/D");
   fChain->Branch("higgsrec_theta",&higgsrec_theta,"higgsrec_theta/D");
   fChain->Branch("higgsrec_E",&higgsrec_E,"higgsrec_E/D");
   fChain->Branch("higgsrec_phi",&higgsrec_phi,"higgsrec_phi/D");
   fChain->Branch("higgsrecem_m",&higgsrecem_m,"higgsrecem_m/D");
   fChain->Branch("higgsrecem_theta",&higgsrecem_theta,"higgsrecem_theta/D");
   fChain->Branch("higgsrecem_E",&higgsrecem_E,"higgsrecem_E/D");
   fChain->Branch("higgsrecem_phi",&higgsrecem_phi,"higgsrecem_phi/D");
   fChain->Branch("higgsHybErec_m",&higgsHybErec_m,"higgsHybErec_m/D");
   fChain->Branch("higgsHybErec_theta",&higgsHybErec_theta,"higgsHybErec_theta/D");
   fChain->Branch("higgsHybErec_E",&higgsHybErec_E,"higgsHybErec_E/D");
   fChain->Branch("higgsHybErec_phi",&higgsHybErec_phi,"higgsHybErec_phi/D");
   fChain->Branch("higgsHybDirRec_m",&higgsHybDirRec_m,"higgsHybDirRec_m/D");
   fChain->Branch("higgsHybDirRec_theta",&higgsHybDirRec_theta,"higgsHybDirRec_theta/D");
   fChain->Branch("higgsHybDirRec_E",&higgsHybDirRec_E,"higgsHybDirRec_E/D");
   fChain->Branch("higgsHybDirRec_phi",&higgsHybDirRec_phi,"higgsHybDirRec_phi/D");
   fChain->Branch("higgstruth_m",&higgstruth_m,"higgstruth_m/D");
   fChain->Branch("higgstruth_theta",&higgstruth_theta,"higgstruth_theta/D");
   fChain->Branch("higgstruth_E",&higgstruth_E,"higgstruth_E/D");
   fChain->Branch("higgstruth_phi",&higgstruth_phi,"higgstruth_phi/D");
   fChain->Branch("higgsscint_m",&higgsscint_m,"higgsscint_m/D");
   fChain->Branch("higgsscint_theta",&higgsscint_theta,"higgsscint_theta/D");
   fChain->Branch("higgsscint_E",&higgsscint_E,"higgsscint_E/D");
   fChain->Branch("higgsscint_phi",&higgsscint_phi,"higgsscint_phi/D");
   fChain->Branch("higgscher_m",&higgscher_m,"higgscher_m/D");
   fChain->Branch("higgscher_theta",&higgscher_theta,"higgscher_theta/D");
   fChain->Branch("higgscher_E",&higgscher_E,"higgscher_E/D");
   fChain->Branch("higgscher_phi",&higgscher_phi,"higgscher_phi/D");
   fChain->Branch("mjjr",&mjjr,"mjjr/D");
   fChain->Branch("mjjt",&mjjt,"mjjt/D");
   fChain->Branch("edep",&edep,"edep/D");
   fChain->Branch("muene_sci",&muene_sci,"muene_sci/D");
   fChain->Branch("muene_che",&muene_che,"muene_che/D");
   fChain->Branch("emcomp1",&emcomp1,"emcomp1/D");
   fChain->Branch("emcomp2",&emcomp2,"emcomp2/D");
   fChain->Branch("etotjr1",&etotjr1,"etotjr1/D");
   fChain->Branch("etotjr2",&etotjr2,"etotjr2/D");
   fChain->Branch("eleak",&eleak,"eleak/D");
   fChain->Branch("eleakn",&eleakn,"eleakn/D");
   fChain->Branch("drmmu",&drmmu,"drmmu/D");
   fChain->Branch("enumu",&enumu,"enumu/D");
   fChain->Branch("mnumu",&mnumu,"mnumu/D");
   fChain->Branch("emu",&emu,"emu/D");
   fChain->Branch("mbos_noc",&mbos_noc,"mbos_noc/D");

   fChain->Branch("j1t_E",&j1t_E,"j1t_E/D");
   fChain->Branch("j1t_pt",&j1t_pt,"j1t_pt/D");
   fChain->Branch("j1t_eta",&j1t_eta,"j1t_eta/D");
   fChain->Branch("j1t_phi",&j1t_phi,"j1t_phi/D");
   fChain->Branch("j1t_m",&j1t_m,"j1t_m/D");
   fChain->Branch("j1t_theta",&j1t_theta,"j1t_theta/D");
   fChain->Branch("j2t_E",&j2t_E,"j2t_E/D");
   fChain->Branch("j2t_pt",&j2t_pt,"j2t_pt/D");
   fChain->Branch("j2t_eta",&j2t_eta,"j2t_eta/D");
   fChain->Branch("j2t_phi",&j2t_phi,"j2t_phi/D");
   fChain->Branch("j2t_m",&j2t_m,"j2t_m/D");
   fChain->Branch("j2t_theta",&j2t_theta,"j2t_theta/D");
   fChain->Branch("j1r_E",&j1r_E,"j1r_E/D");
   fChain->Branch("j1r_pt",&j1r_pt,"j1r_pt/D");
   fChain->Branch("j1r_eta",&j1r_eta,"j1r_eta/D");
   fChain->Branch("j1r_phi",&j1r_phi,"j1r_phi/D");
   fChain->Branch("j1r_m",&j1r_m,"j1r_m/D");
   fChain->Branch("j1r_theta",&j1r_theta,"j1r_theta/D");
   fChain->Branch("j2r_E",&j2r_E,"j2r_E/D");
   fChain->Branch("j2r_pt",&j2r_pt,"j2r_pt/D");
   fChain->Branch("j2r_eta",&j2r_eta,"j2r_eta/D");
   fChain->Branch("j2r_phi",&j2r_phi,"j2r_phi/D");
   fChain->Branch("j2r_m",&j2r_m,"j2r_m/D");
   fChain->Branch("j2r_theta",&j2r_theta,"j2r_theta/D");
   fChain->Branch("j1s_E",&j1s_E,"j1s_E/D");
   fChain->Branch("j1s_pt",&j1s_pt,"j1s_pt/D");
   fChain->Branch("j1s_eta",&j1s_eta,"j1s_eta/D");
   fChain->Branch("j1s_phi",&j1s_phi,"j1s_phi/D");
   fChain->Branch("j1s_m",&j1s_m,"j1s_m/D");
   fChain->Branch("j1s_theta",&j1s_theta,"j1s_theta/D");
   fChain->Branch("j2s_E",&j2s_E,"j2s_E/D");
   fChain->Branch("j2s_pt",&j2s_pt,"j2s_pt/D");
   fChain->Branch("j2s_eta",&j2s_eta,"j2s_eta/D");
   fChain->Branch("j2s_phi",&j2s_phi,"j2s_phi/D");
   fChain->Branch("j2s_m",&j2s_m,"j2s_m/D");
   fChain->Branch("j2s_theta",&j2s_theta,"j2s_theta/D");
   fChain->Branch("j1c_E",&j1c_E,"j1c_E/D");
   fChain->Branch("j1c_pt",&j1c_pt,"j1c_pt/D");
   fChain->Branch("j1c_eta",&j1c_eta,"j1c_eta/D");
   fChain->Branch("j1c_phi",&j1c_phi,"j1c_phi/D");
   fChain->Branch("j1c_m",&j1c_m,"j1c_m/D");
   fChain->Branch("j1c_theta",&j1c_theta,"j1c_theta/D");
   fChain->Branch("j2c_E",&j2c_E,"j2c_E/D");
   fChain->Branch("j2c_pt",&j2c_pt,"j2c_pt/D");
   fChain->Branch("j2c_eta",&j2c_eta,"j2c_eta/D");
   fChain->Branch("j2c_phi",&j2c_phi,"j2c_phi/D");
   fChain->Branch("j2c_m",&j2c_m,"j2c_m/D");
   fChain->Branch("j2c_theta",&j2c_theta,"j2c_theta/D");
   fChain->Branch("closestT_DR_gam1",&closestT_DR_gam1,"closestT_DR_gam1/D");
   fChain->Branch("closestT_DR_gam2",&closestT_DR_gam2,"closestT_DR_gam2/D");

}

void MyTree::Reset()
{
   nmuon =-1;
   nneu =-1;
   mjjr=-1.;
   mjjt=-1.;
   edep=-1.;
   muene_sci=-1.;
   muene_che=-1.;
   emcomp1=-1.;
   emcomp2=-1.;
   etotjr1=-1.;
   etotjr2=-1.;
   drmmu=-1.;

   j1t_E=-1.;
   j1t_pt=-1.;
   j1t_eta=-10.;
   j1t_phi=-5.;
   j1t_m=-1.;
   j1t_theta=-7.;
   j2t_E=-1.;
   j2t_pt=-1.;
   j2t_eta=-10.;
   j2t_phi=-5.;
   j2t_m=-1.;
   j2t_theta=-7.;
   j1r_E=-1.;
   j1r_pt=-1.;
   j1r_eta=-10.;
   j1r_phi=-5.;
   j1r_m=-1.;
   j1r_theta=-7.;
   j2r_E=-1.;
   j2r_pt=-1.;
   j2r_eta=-10.;
   j2r_phi=-5.;
   j2r_m=-1.;
   j2r_theta=-7.;
   j1s_E=-1.;
   j1s_pt=-1.;
   j1s_eta=-10.;
   j1s_phi=-5.;
   j1s_m=-1.;
   j1s_theta=-7.;
   j2s_E=-1.;
   j2s_pt=-1.;
   j2s_eta=-10.;
   j2s_phi=-5.;
   j2s_m=-1.;
   j2s_theta=-7.;
   j1c_E=-1.;
   j1c_pt=-1.;
   j1c_eta=-10.;
   j1c_phi=-5.;
   j1c_m=-1.;
   j1c_theta=-7.;
   j2c_E=-1.;
   j2c_pt=-1.;
   j2c_eta=-10.;
   j2c_phi=-5.;
   j2c_m=-1.;
   j2c_theta=-7.;

   higgsrec_m = 0;
   higgsrec_theta = 0;
   higgsrec_E = 0;
   higgsrec_phi = 0;
   higgsrecem_m = 0;
   higgsrecem_theta = 0;
   higgsrecem_E = 0;
   higgsrecem_phi = 0;
   higgsHybErec_m = 0;
   higgsHybErec_theta = 0;
   higgsHybErec_E = 0;
   higgsHybErec_phi = 0;
   higgsHybDirRec_m = 0;
   higgsHybDirRec_theta = 0;
   higgsHybDirRec_E = 0;
   higgsHybDirRec_phi = 0;
   higgstruth_m = 0;
   higgstruth_theta = 0;
   higgstruth_E = 0;
   higgstruth_phi = 0;
   higgsscint_m = 0;
   higgsscint_theta = 0;
   higgsscint_E = 0;
   higgsscint_phi = 0;
   higgscher_m = 0;
   higgscher_theta = 0;
   higgscher_E = 0;
   higgscher_phi = 0;

}
