/***********************************************************************
	This is the main program to process the CBNT Root Ntuple from
Athena with SUSYtup. See SUSYtup.h for more information.
	For version 7.0.0++ using name susy for Ntuple.
***********************************************************************/
// ROOT includes 

#include <TROOT.h> 
#include "TTree.h" 
#include "TBranch.h" 
#include "TProfile.h"
#include <TFile.h>
#include "TLorentzVector.h"
#include "TMath.h"

// fastjet includes 

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

// Framework includes

#include "MyTree.h"

// std library includes 

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <algorithm>

using namespace std;

TFile* ftree;
MyTree bonsaiTree;

vector<fastjet::PseudoJet> inputparticles_tru;
vector<fastjet::PseudoJet> inputparticles_scin;
vector<fastjet::PseudoJet> inputparticles_cher;
vector<fastjet::PseudoJet> jetexc;
vector<fastjet::PseudoJet> jet_scin;
vector<fastjet::PseudoJet> jet_cher_t;
vector<fastjet::PseudoJet> jet_cher;
vector<fastjet::PseudoJet> jet_rec;
vector<fastjet::PseudoJet> jet_recem; // reco jet taking simply the average of the two signals
vector<fastjet::PseudoJet> jet_tru;
vector<fastjet::PseudoJet> jet_truem;
vector<TLorentzVector> gamvec;
vector<TLorentzVector> nuvec;
vector<TLorentzVector> otherparticle;
vector <double> emcomp;
vector <double> etotjr;

vector<double> Calib_VectorScinR;
vector<double> Calib_VectorScinL;
vector<double> Calib_VectorCherR;
vector<double> Calib_VectorCherL;	

double GeV=1000.;
double pi=TMath::Pi();
double threshold=0.01;
double etalim=5; //5.0

//////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) { 
  void recalibrate(std::vector<fastjet::PseudoJet> * jet, TProfile * h_cal);
  fastjet::PseudoJet findClosestTower(fastjet::PseudoJet v);  
  std::vector<unsigned int> cutflow(10,0.);  
  vector<double> calibscin(std::vector<double> vectorscin);
  vector<double> calibcher(std::vector<double> vectorcher);
  tuple<double, double, double> maptower(int index, string side);
  fastjet::PseudoJet makeHybrid(fastjet::PseudoJet mag, fastjet::PseudoJet direction);
  fastjet::PseudoJet matchjet(fastjet::PseudoJet jet_in, vector<fastjet::PseudoJet> testvec);
  fastjet::PseudoJet mergejet(fastjet::PseudoJet jet_scin, fastjet::PseudoJet jet_cher);
  fastjet::PseudoJet mergejetem(fastjet::PseudoJet jet_scin, fastjet::PseudoJet jet_cher, bool doCalib); 
  double energyShare(fastjet::PseudoJet jet);
  cout << "-------------------------------------------" << std::endl;
  cout << " N arguments " << argc << std::endl;
  if(argc<2){
    cout << "Please give output.root input.root" << endl;
    return 0;
  }


  bool doRecal = true;

  TProfile * h_calib_scint = 0;
  TProfile * h_calib_cher = 0;
  TFile * fCal = 0;
  if (doRecal){
    fCal = TFile::Open("calfile_DR_hgamgam.root");
    if (!fCal){
      std::cerr << "Cannot find calibration file" << std::endl;
    } else {
      h_calib_scint = (TProfile *) fCal->Get("h_calib_scint1");
      h_calib_cher =  (TProfile *) fCal->Get("h_calib_cher1");
    }
  }
    
  
  // Parsing of the output file
  
  // Output tree file
  string histName=argv[1];
  std::cout << " Output file name: " << std::endl;
  std::cout << "      " << histName << std::endl;
  
  // Open input file
  std::string fn = argv[2];
  std::cout << " Read file name: " << std::endl;
  std::cout << "      " << fn << std::endl;
  
  // Making the assumption on the name of the output file. Maybe the two trees should at some point be written in the same output file? 
  
  std::string filetru = fn+"_truth.root";
  std::string filesim = fn+".root";
  std::cout<< "Reading file " <<filetru <<std::endl;
  TFile* f = new TFile( filetru.c_str() );
  std::cout<< "Reading file " <<filesim <<std::endl;
  TFile* f1 = new TFile( filesim.c_str() );
  
  int pos_st=histName.rfind("/");
  string stma=histName.substr(pos_st+1);
  string newfile="bonsai."+stma;
  cout << " bonsai file " << newfile << endl;
  ftree = new TFile(newfile.c_str(), "RECREATE");
  bonsaiTree.Init();
  
#include "truthdec.h"
#include "B4dec.h"
  
  //
  TTree* tree1 = (TTree*)f->Get("truth");
  TTree* tree2 = (TTree*)f1->Get("B4");
  tree1->AddFriend(tree2);
#include "truthset.h"
#include "B4set.h"
  tree1->GetEntry(0);
  //
  
  
  if (tree1== 0) return 1;
  Int_t nentries = Int_t(tree1->GetEntries());
  Int_t nbytes= 0, nb = 0;
  // Loop on the number of events
  
  cout << " Number of events " << nentries << endl;
  for (Int_t jentry=0; jentry<nentries;jentry++) {
    if (jentry % 100 == 0) std::cout << "Processed " << jentry << " entries" << std::endl;
    nb = tree1->GetEntry(jentry);   nbytes += nb;
    bonsaiTree.Reset();
    inputparticles_tru.clear();
    gamvec.clear();
    nuvec.clear();
    int ngam=0;
    int nneu=0;
    int nmun=0;
    double gamene_sci=0.;
    double gamene_che=0.;
    double etott=0;
    
    // Now looping over the truth particles 
    // Now the event reconstruction is finished. 
    
    // Check how many events we have at different stages of teh selection
    
    ++cutflow[0];
    
    for(uint itru=0;itru<mcs_n;itru++){
      int partid = mcs_pdgId->at(itru);
      double parteta = mcs_eta->at(itru);
      etott+=mcs_E->at(itru);
      
      // Prepare them to be clustered by fastjet if they are not neutrinos or neutralinos, or muons
      
      if(abs(partid) != 13 &&  
         abs(partid) !=12  && abs(partid) != 14 && abs(partid) != 16 &&
         abs(partid) != 1000022){
        if(abs(parteta)<etalim){
          TLorentzVector trup;
          trup.SetPtEtaPhiM(mcs_pt->at(itru), mcs_eta->at(itru), mcs_phi->at(itru),
                            mcs_m->at(itru));
          fastjet::PseudoJet fj(trup.Px(), trup.Py(), trup.Pz(), trup.E());
          fj.set_user_index(itru);
          inputparticles_tru.push_back(fj);
        }
      }
      
      if(abs(partid) ==12  ||  abs(partid) == 14 || abs(partid) == 16){
        TLorentzVector nu_p;
        nu_p.SetPtEtaPhiM(mcs_pt->at(itru), mcs_eta->at(itru), mcs_phi->at(itru),mcs_m->at(itru));
        nuvec.push_back(nu_p);
	nneu++;
      }
      
      if(abs(partid) == 22){
        TLorentzVector gam;
        gam.SetPtEtaPhiM(mcs_pt->at(itru), mcs_eta->at(itru), mcs_phi->at(itru),
			 mcs_m->at(itru));
        gamvec.push_back(gam);
        ngam++;
      } 
    } // loop on truth particles    
    //    cout << etott << endl;
    jetexc.clear();
    
    if (mcs_n != 4) continue; 
    ++cutflow[1];

    if (ngam != 2) continue;
    ++cutflow[2];
    
    // Define what jet should be reconstructed
    
    fastjet::JetDefinition jet_def(fastjet::ee_genkt_algorithm, 2.*pi, 1.);
    fastjet::ClusterSequence clust_seq(gamvec, jet_def); 
    
    // This will be a jet sequence containing the "truth" jets 
    
    jetexc = clust_seq.exclusive_jets(int(2));;
    
    if (jetexc.size() != 2){
      std::cout << "Why is this different from 2?" << std::endl;
      continue;
    }
    
    // Check about the DR separation of the two photons. Reject events where the two photons are not well separated.
    
    if (jetexc[0].delta_R(jetexc[1]) < 1) continue;
    ++cutflow[3];
    
    //  now the rec part
    //
    Calib_VectorScinR.clear();
    Calib_VectorScinL.clear();
    Calib_VectorCherR.clear();
    Calib_VectorCherL.clear();
    // 
    Calib_VectorScinR = calibscin(*VectorSignalsR);
    Calib_VectorScinL = calibscin(*VectorSignalsL);
    Calib_VectorCherR = calibcher(*VectorSignalsCherR);
    Calib_VectorCherL = calibcher(*VectorSignalsCherL);
    
    double energy=0;
    
    for(uint i=0; i<Calib_VectorScinR.size(); i++) {
      energy+=Calib_VectorScinR.at(i)+Calib_VectorScinL.at(i);
    }       
    
    if(energy>0){
      inputparticles_scin.clear();
      inputparticles_cher.clear();
      TLorentzVector scin_bos(0.,0.,0.,0.);
      TLorentzVector cher_bos(0.,0.,0.,0.);
      // right side
      for(int towerindex=1; towerindex<=75*36; towerindex++) {
        auto thphieta=maptower(towerindex, "r");
        double theta=get<0>(thphieta);
        double phi=get<1>(thphieta);
        double eta=get<2>(thphieta);
        double energy_scin = Calib_VectorScinR[towerindex];
        double pt_scin = energy_scin*sin(theta*pi/180.);
        double energy_cher = Calib_VectorCherR[towerindex];
        double pt_cher = energy_cher*sin(theta*pi/180.);
        TLorentzVector towerscin;
        towerscin.SetPtEtaPhiM(pt_scin, eta, phi*pi/180., 0.);
        TLorentzVector towercher;
        towercher.SetPtEtaPhiM(pt_cher, eta, phi*pi/180., 0.);
	if (TMath::Abs(towerscin.E()) > 0.01) inputparticles_scin.push_back(fastjet::PseudoJet(towerscin.Px(),towerscin.Py(), towerscin.Pz(), towerscin.E()));
	if (TMath::Abs(towercher.E()) > 0.01) inputparticles_cher.push_back(fastjet::PseudoJet(towercher.Px(),towercher.Py(), towercher.Pz(), towercher.E()));
      } 
      
      // left sid3
      for(int towerindex=1; towerindex<=75*36; towerindex++) {
        auto thphieta=maptower(towerindex, "l");
        double theta=get<0>(thphieta);
        double phi=get<1>(thphieta);
        double eta=get<2>(thphieta);
	
	//        cout << towerindex << " theta " << theta << " phi " << phi << " eta " << eta << endl;
        double energy_scin = Calib_VectorScinL[towerindex];
        double pt_scin = energy_scin*sin(theta*pi/180.);
	
        double energy_cher = Calib_VectorCherL[towerindex];
        double pt_cher = energy_cher*sin(theta*pi/180.);
	
        TLorentzVector towerscin;
        towerscin.SetPtEtaPhiM(pt_scin, eta, phi*pi/180., 0.);
        TLorentzVector towercher;
        towercher.SetPtEtaPhiM(pt_cher, eta, phi*pi/180., 0.);
	if (TMath::Abs(towerscin.E()) > 0.01) inputparticles_scin.push_back(fastjet::PseudoJet(towerscin.Px(),towerscin.Py(), towerscin.Pz(), towerscin.E()));
	if (TMath::Abs(towercher.E()) > 0.01) inputparticles_cher.push_back(fastjet::PseudoJet(towercher.Px(),towercher.Py(), towercher.Pz(), towercher.E()));
      }
      
      //
      /*      fastjet::JetDefinition jet_defs(fastjet::ee_genkt_algorithm, 2.*pi, 1.);
      fastjet::ClusterSequence clust_seq_scin(inputparticles_scin, jet_defs); 
      fastjet::ClusterSequence clust_seq_cher(inputparticles_cher, jet_defs);
      */


      fastjet::JetDefinition jet_defs(fastjet::ee_genkt_algorithm, 2.*pi, 1.);
      fastjet::ClusterSequence clust_seq_scin(inputparticles_scin, jet_defs); 
      fastjet::ClusterSequence clust_seq_cher(inputparticles_cher, jet_defs);
      
      //  clear vectors of jets
      jet_scin.clear();
      jet_cher.clear();
      jet_cher_t.clear();
      jet_rec.clear();
      jet_recem.clear();
      jet_tru.clear();
      jet_truem.clear();
      //  create vector of jet_scin
      jet_scin = clust_seq_scin.exclusive_jets(int(2));;
      
      
      //   create temp vector of jet_cher 
      jet_cher_t = clust_seq_cher.exclusive_jets(int(2));;
      
      /*      std::cout << "New event" << std::endl;
      std::cout << "NJets scint " << jet_scin.size() << std::endl;
      for (unsigned int i = 0; i < jet_scin.size(); ++i){
	std::cout << "jet_scin " << i << " Pt, Eta, Phi " << jet_scin[i].pt() << "   " << jet_scin[i].eta() << "   " << jet_scin[i].phi() << std::endl;
      }
      std::cout << "NJets cher " << jet_cher_t.size() << std::endl;
      for (unsigned int i = 0; i < jet_cher_t.size(); ++i){
	std::cout << "jet_cher_t " << i << " Pt, Eta, Phi " << jet_cher_t[i].pt() << "   " << jet_cher_t[i].eta() << "   " << jet_cher_t[i].phi() << std::endl;
	}*/

      
      //   align jet_cher and jet_scin vector
      for(uint jn=0; jn<jet_scin.size();jn++) {
        jet_cher.push_back(matchjet(jet_scin[jn], jet_cher_t)); 
      }

      // Recalibrate the jets if requested

      if (doRecal){
	recalibrate(&jet_scin,h_calib_scint);
	recalibrate(&jet_cher,h_calib_cher);
      }

      //   combine aligned scin and cher into recem
      for(uint jn=0; jn<jet_scin.size();jn++) {
        jet_recem.push_back(mergejetem(jet_scin[jn],jet_cher[jn],doRecal));
      }

      for(uint jn=0; jn<jet_recem.size();jn++) {
        jet_truem.push_back(matchjet(jet_recem[jn], jetexc)); 
      }

      //   combine aligned scin and cher into rec
      for(uint jn=0; jn<jet_scin.size();jn++) {
        jet_rec.push_back(mergejet(jet_scin[jn],jet_cher[jn]));
      }
      
      for(uint jn=0; jn<jet_rec.size();jn++) {
        jet_tru.push_back(matchjet(jet_rec[jn], jetexc)); 
      }
      
      if (jet_scin.size() != 2) continue;
      ++cutflow[4];
      if (jet_cher.size() != 2) continue;
      ++cutflow[5];
      if (jet_rec.size() != 2) continue;
      ++cutflow[6];
      if (jet_recem.size() != 2) continue;
      ++cutflow[7];
      if (jet_recem[0].E() < 10 || jet_recem[1].E() < 10 || jet_tru[0].E() < 10 || jet_tru[1].E() < 10) continue;
      ++cutflow[8];
    

      //
      //   muiso for W
      //
      double drminmu=9999.;
      if(gamvec.size()>0) {	      
	for(uint jt=0; jt<jet_tru.size(); jt++) {
	  TLorentzVector jj;
	  jj.SetPtEtaPhiM(jet_tru[jt].pt(), jet_tru[jt].eta(), jet_tru[jt].phi(),
			  jet_tru[jt].m());
	  double deltaR=abs(jj.DeltaR(gamvec[0]));
	  if(deltaR<drminmu)drminmu=deltaR;
	}
      }

      double emu=0.;
      double enumu=0; 
      double mnumu=0;
      if(ngam==1)emu=gamvec[0].E();
      if(ngam==1 && nmun==1) {
	enumu=gamvec[0].E()+nuvec[0].E();
	mnumu=(gamvec[0]+nuvec[0]).M();
      }
      //
      //    save in ntuple
      //
      
      fastjet::PseudoJet jetrec(0.,0.,0.,0.);
      fastjet::PseudoJet jettruth(0.,0.,0.,0.);
      fastjet::PseudoJet jetscint(0.,0.,0.,0.);
      fastjet::PseudoJet jetcher(0.,0.,0.,0.);
      fastjet::PseudoJet jetrecem(0.,0.,0.,0.);
      fastjet::PseudoJet jetrecem_recal(0.,0.,0.,0.);
      fastjet::PseudoJet jethybErec(0.,0.,0.,0.);
      fastjet::PseudoJet jethybDirRec(0.,0.,0.,0.);
      if(jet_rec.size()==2 && jet_tru.size()==2) {
        jetrec=jet_rec[0]+jet_rec[1];
	jetrecem=jet_recem[0] + jet_recem[1];
        jettruth=jet_tru[0]+jet_tru[1];
	jetscint = jet_scin[0] + jet_scin[1];
	jetcher = jet_cher[0] + jet_cher[1];
	jethybErec = makeHybrid(jet_recem[0],jet_tru[0]) + makeHybrid(jet_recem[1],jet_tru[1]);
	jethybDirRec = makeHybrid(jet_tru[0],jet_recem[0]) + makeHybrid(jet_tru[1],jet_recem[1]);
      }

      double eshare_scin0 = energyShare(jet_scin[0]);      
      double eshare_scin1 = energyShare(jet_scin[1]);      
      double eshare_cher0 = energyShare(jet_cher[0]);      
      double eshare_cher1 = energyShare(jet_cher[1]);      


      bonsaiTree.nmuon = ngam;
      bonsaiTree.nneu = nneu;
      bonsaiTree.higgsrec_m= jetrec.m();
      bonsaiTree.higgsrec_theta= jetrec.theta();
      bonsaiTree.higgsrec_phi= jetrec.phi();
      bonsaiTree.higgsrec_E= jetrec.E();
      bonsaiTree.higgsrecem_m= jetrecem.m();
      bonsaiTree.higgsrecem_theta= jetrecem.theta();
      bonsaiTree.higgsrecem_phi= jetrecem.phi();
      bonsaiTree.higgsrecem_E= jetrecem.E();
      bonsaiTree.higgsHybErec_m= jethybErec.m();
      bonsaiTree.higgsHybErec_theta= jethybErec.theta();
      bonsaiTree.higgsHybErec_phi= jethybErec.phi();
      bonsaiTree.higgsHybErec_E= jethybErec.E();
      bonsaiTree.higgsHybDirRec_m= jethybDirRec.m();
      bonsaiTree.higgsHybDirRec_theta= jethybDirRec.theta();
      bonsaiTree.higgsHybDirRec_phi= jethybDirRec.phi();
      bonsaiTree.higgsHybDirRec_E= jethybDirRec.E();
      bonsaiTree.higgstruth_theta= jettruth.theta();
      bonsaiTree.higgstruth_phi= jettruth.phi();
      bonsaiTree.higgstruth_E= jettruth.E();
      bonsaiTree.higgstruth_m= jettruth.m();
      bonsaiTree.higgsscint_theta= jetscint.theta();
      bonsaiTree.higgsscint_phi= jetscint.phi();
      bonsaiTree.higgsscint_E= jetscint.E();
      bonsaiTree.higgsscint_m= jetscint.m();
      bonsaiTree.higgscher_theta= jetcher.theta();
      bonsaiTree.higgscher_phi= jetcher.phi();
      bonsaiTree.higgscher_E= jetcher.E();
      bonsaiTree.higgscher_m= jetcher.m();
      bonsaiTree.edep= EnergyTot/1000.;
      bonsaiTree.muene_sci=gamene_sci;
      bonsaiTree.muene_che=gamene_che;
      bonsaiTree.emcomp1=(emcomp.size() > 0) ? emcomp[0] : -10000;
      bonsaiTree.emcomp2=(emcomp.size() > 1) ? emcomp[1] : -10000;
      bonsaiTree.etotjr1=(etotjr.size() > 0) ? etotjr[0] : -10000;
      bonsaiTree.etotjr2=(etotjr.size() > 1) ? etotjr[1] : -10000;
      bonsaiTree.eleak=leakage;
      bonsaiTree.eleakn=neutrinoleakage;
      bonsaiTree.mbos_noc=-1;
      bonsaiTree.drmmu=drminmu;
      bonsaiTree.enumu=enumu;
      bonsaiTree.mnumu=mnumu;
      bonsaiTree.emu=emu;
      
      //
      
      bonsaiTree.j1t_E=(jet_tru.size() > 0) ? jet_tru[0].E() : -10000;
      bonsaiTree.j1t_pt=(jet_tru.size() > 0) ? jet_tru[0].pt() : -10000;
      bonsaiTree.j1t_eta=(jet_tru.size() > 0) ? jet_tru[0].eta() : -10000;
      bonsaiTree.j1t_phi=(jet_tru.size() > 0) ? jet_tru[0].phi() : -10000;
      bonsaiTree.j1t_m=(jet_tru.size() > 0) ? jet_tru[0].m() : -10000;
      bonsaiTree.j1t_theta=(jet_tru.size() > 0) ? jet_tru[0].theta() : -10000;
      bonsaiTree.j2t_E=(jet_tru.size() > 1) ? jet_tru[1].E() : -10000;
      bonsaiTree.j2t_pt=(jet_tru.size() > 1) ? jet_tru[1].pt() : -10000;
      bonsaiTree.j2t_eta=(jet_tru.size() > 1) ? jet_tru[1].eta() : -10000;
      bonsaiTree.j2t_phi=(jet_tru.size() > 1) ? jet_tru[1].phi() : -10000;
      bonsaiTree.j2t_m=(jet_tru.size() > 1) ? jet_tru[1].m() : -10000;
      bonsaiTree.j2t_theta=(jet_tru.size() > 1) ? jet_tru[1].theta() : -10000;
      bonsaiTree.j1r_E=(jet_rec.size() > 0) ? jet_rec[0].E() : -10000;
      bonsaiTree.j1r_pt=(jet_rec.size() > 0) ? jet_rec[0].pt() : -10000;
      bonsaiTree.j1r_eta=(jet_rec.size() > 0) ? jet_rec[0].eta() : -10000;
      bonsaiTree.j1r_phi=(jet_rec.size() > 0) ? jet_rec[0].phi() : -10000;
      bonsaiTree.j1r_m=(jet_rec.size() > 0) ? jet_rec[0].m() : -10000;
      bonsaiTree.j1r_theta=(jet_rec.size() > 0) ? jet_rec[0].theta() : -10000;
      bonsaiTree.j2r_E=(jet_rec.size() > 1) ? jet_rec[1].E() : -10000;
      bonsaiTree.j2r_pt=(jet_rec.size() > 1) ? jet_rec[1].pt() : -10000;
      bonsaiTree.j2r_eta=(jet_rec.size() > 1) ? jet_rec[1].eta() : -10000;
      bonsaiTree.j2r_phi=(jet_rec.size() > 1) ? jet_rec[1].phi() : -10000;
      bonsaiTree.j2r_m=(jet_rec.size() > 1) ? jet_rec[1].m() : -10000;
      bonsaiTree.j2r_theta=(jet_rec.size() > 1) ? jet_rec[1].theta() : -10000;
      bonsaiTree.j1s_E=(jet_scin.size() > 0) ? jet_scin[0].E() : -10000;
      bonsaiTree.j1s_pt=(jet_scin.size() > 0) ? jet_scin[0].pt() : -10000;
      bonsaiTree.j1s_eta=(jet_scin.size() > 0) ? jet_scin[0].eta() : -10000;
      bonsaiTree.j1s_phi=(jet_scin.size() > 0) ? jet_scin[0].phi() : -10000;
      bonsaiTree.j1s_m=(jet_scin.size() > 0) ? jet_scin[0].m() : -10000;
      bonsaiTree.j1s_theta=(jet_scin.size() > 0) ? jet_scin[0].theta() : -10000;
      bonsaiTree.j2s_E=(jet_scin.size() > 1) ? jet_scin[1].E() : -10000;
      bonsaiTree.j2s_pt=(jet_scin.size() > 1) ? jet_scin[1].pt() : -10000;
      bonsaiTree.j2s_eta=(jet_scin.size() > 1) ? jet_scin[1].eta() : -10000;
      bonsaiTree.j2s_phi=(jet_scin.size() > 1) ? jet_scin[1].phi() : -10000;
      bonsaiTree.j2s_m=(jet_scin.size() > 1) ? jet_scin[1].m() : -10000;
      bonsaiTree.j2s_theta=(jet_scin.size() > 1) ? jet_scin[1].theta() : -10000;
      bonsaiTree.j1c_E=(jet_cher.size() > 0) ? jet_cher[0].E() : -10000;
      bonsaiTree.j1c_pt=(jet_cher.size() > 0) ? jet_cher[0].pt() : -10000;
      bonsaiTree.j1c_eta=(jet_cher.size() > 0) ? jet_cher[0].eta() : -10000;
      bonsaiTree.j1c_phi=(jet_cher.size() > 0) ? jet_cher[0].phi() : -10000;
      bonsaiTree.j1c_m=(jet_cher.size() > 0) ? jet_cher[0].m() : -10000;
      bonsaiTree.j1c_theta=(jet_cher.size() > 0) ? jet_cher[0].theta() : -10000;
      bonsaiTree.j2c_E=(jet_cher.size() > 1) ? jet_cher[1].E() : -10000;
      bonsaiTree.j2c_pt=(jet_cher.size() > 1) ? jet_cher[1].pt() : -10000;
      bonsaiTree.j2c_eta=(jet_cher.size() > 1) ? jet_cher[1].eta() : -10000;
      bonsaiTree.j2c_phi=(jet_cher.size() > 1) ? jet_cher[1].phi() : -10000;
      bonsaiTree.j2c_m=(jet_cher.size() > 1) ? jet_cher[1].m() : -10000;
      bonsaiTree.j2c_theta=(jet_cher.size() > 1) ? jet_cher[1].theta() : -10000;
      bonsaiTree.j1s_eshare = eshare_scin0;
      bonsaiTree.j2s_eshare = eshare_scin1;
      bonsaiTree.j1c_eshare = eshare_cher0;
      bonsaiTree.j2c_eshare = eshare_cher1;
      //      bonsaiTree.closestT_DR_gam1 = (jet_tru.size() > 0) ?  jet_tru[0].delta_R(findClosestTower(jet_tru[0])) : -10000;
      //      bonsaiTree.closestT_DR_gam2 = (jet_tru.size() > 1) ?  jet_tru[1].delta_R(findClosestTower(jet_tru[1])) : -10000;
      bonsaiTree.closestT_DR_gam1 = 0;
      bonsaiTree.closestT_DR_gam2 = 0;
          
    }// energy>0          
    
// fill output tree
    bonsaiTree.Fill();
  } // loop on events
  
  for (unsigned int i = 0; i < cutflow.size(); ++i){
    std::cout << "cutflow[" << i << "] = " << cutflow[i] << std::endl;
  }

  ftree->cd();
  bonsaiTree.Write();
  delete f;
  delete f1;

}





//
std::tuple<double, double, double> maptower(int index, std::string side){
//Function to return tower angles (theta and phi) given index
  int NbOfBarrel=40;
  int NbOfEndcap=35;
  int NZrot=36;
  int TotTower=NbOfBarrel+NbOfEndcap;
  index = index-1;
  int sliceindex = index/TotTower;
  int towerindex = index-(sliceindex*TotTower);
  double deltatheta = 45./(NbOfBarrel);
//get theta
  double theta = towerindex*deltatheta+deltatheta/2.;
//  cout << " thetap " << theta << endl;
//get phi
  double phi_unit = 360./NZrot;
  double phi = (sliceindex)*phi_unit;
  
  if (side == "r"){
//     cout << " thetai " << theta+90. << " phii " << phi << " etai " <<  -log(tan(((90.-theta)*pi/180./2.))) << endl;
    return std::make_tuple(theta+90., phi, -log(tan(((90.-theta)*pi/180./2.))));
  }
  if (side == "l"){
    return std::make_tuple(90.-theta, phi, log(tan(((90.-theta)*pi/180./2.))));
  }
  return std::make_tuple(0.,0.,0.);
}
fastjet::PseudoJet mergejet(fastjet::PseudoJet jet_scin, fastjet::PseudoJet jet_cher) {

  double c=0.34;

  double jetPx = (jet_scin.px()-c*jet_cher.px())/(1-c);
  double jetPy = (jet_scin.py()-c*jet_cher.py())/(1-c);
  double jetPz = (jet_scin.pz()-c*jet_cher.pz())/(1-c);
  double jetE = (jet_scin.e()-c*jet_cher.e())/(1.-c);
  return fastjet::PseudoJet(jetPx, jetPy, jetPz, jetE);
}

fastjet::PseudoJet mergejetem(fastjet::PseudoJet jet_scin, fastjet::PseudoJet jet_cher, bool doCalib) {

  double cherConst = 1.;
  double scinConst = 1.;

  if (doCalib){
    cherConst = 1;
    scinConst = 1;
  }
  
  double jetPx = (scinConst*jet_scin.px()+ cherConst*jet_cher.px())/2.;
  double jetPy = (scinConst*jet_scin.py()+ cherConst*jet_cher.py())/2.;
  double jetPz = (scinConst*jet_scin.pz()+ cherConst*jet_cher.pz())/2.;
  double jetE = (scinConst*jet_scin.e() + cherConst*jet_cher.e())/2.;
  return fastjet::PseudoJet(jetPx, jetPy, jetPz, jetE);
}

fastjet::PseudoJet matchjet(fastjet::PseudoJet jet_in, vector<fastjet::PseudoJet> testvec) {

  int imin=-1;
  double deltarmin=99999.;
  for(uint i=0; i<testvec.size(); i++){
    double deltar=jet_in.delta_R(testvec.at(i));
    if(deltar<deltarmin) {
      deltarmin=deltar;
      imin=i;
    }
  }
  if(imin != -1) return testvec.at(imin);
  else
  return fastjet::PseudoJet(0., 0., 0., 0.);
}

fastjet::PseudoJet makeHybrid(fastjet::PseudoJet mag, fastjet::PseudoJet direction)
{
  double magnitude = mag.E();
  double dirMag = TMath::Sqrt(TMath::Power(direction.px(),2) +
			      TMath::Power(direction.py(),2) +
			      TMath::Power(direction.pz(),2));

  if (TMath::Abs(dirMag) < 0.000001) {
    std::cout << "Warning null magnitude, returning 0" << std::endl;
    return fastjet::PseudoJet(0., 0., 0., 0.);
  }
			      
  return fastjet::PseudoJet(magnitude * direction.px()/dirMag, magnitude * direction.py()/dirMag, magnitude * direction.pz()/dirMag, magnitude);
}

fastjet::PseudoJet findClosestTower(fastjet::PseudoJet v)
{
  fastjet::PseudoJet retval;
  fastjet::PseudoJet tempval;
  TLorentzVector tlv;
  double DR = 100000;
  for(int towerindex=1; towerindex<=75*36; towerindex++) {                  
    auto thphieta=maptower(towerindex, "r");                            
    double theta=get<0>(thphieta);                                          
    double phi=get<1>(thphieta);                                            
    double eta=get<2>(thphieta);    
    tlv.SetPtEtaPhiM(1., eta, phi*pi/180., 0.);
    tempval.reset(tlv.Px(),tlv.Py(),tlv.Pz(),tlv.E());
    if (tempval.delta_R(v) < DR){
      retval.reset(tempval);
      DR = tempval.delta_R(v);
    }
    thphieta=maptower(towerindex, "l");                            
    theta=get<0>(thphieta);                                          
    phi=get<1>(thphieta);                                            
    eta=get<2>(thphieta);    
    tlv.SetPtEtaPhiM(1., eta, phi*pi/180., 0.);
    tempval.reset(tlv.Px(),tlv.Py(),tlv.Pz(),tlv.E());
    if (tempval.delta_R(v) < DR){
      retval.reset(tempval);
      DR = tempval.delta_R(v);
    }
  }

  return retval;
}

double energyShare(fastjet::PseudoJet jet)
{
  double retval = -1.;

  std::vector<fastjet::PseudoJet> j_const = jet.constituents();
  std::vector<double> energies;
  energies.reserve(j_const.size());
  for (auto i : j_const){
    energies.push_back(i.E());
  }

  std::sort(energies.begin(),energies.end(),std::greater<double>());
  
  if (energies.size() < 2) retval = 0.;
  else retval = energies[1]/energies[0];

  return retval;
}

/************************ original calibration constants ************

std::vector<double> calibscin(std::vector<double> vectorscin){
	std::vector<double> s_cont;
	s_cont = {391.9964784225476, 392.65450934717455, 391.9130470386118, 390.6050630558208, 389.1126112802391, 387.95945028322916, 387.4823598282578, 388.1313680310084, 389.1887061971629, 389.39650703390834, 388.1539715075822, 388.50402528038387, 388.8469797524201, 389.5104813700643, 389.9263913577828, 389.03950145952115, 388.99508641977286, 388.772944352741, 389.1309191354679, 389.2350876534651, 389.13791528990845, 388.8450419930907, 389.1250633143002, 388.9702856267768, 388.89437597892635, 389.54155385650677, 389.6077230144546, 389.83471955121126, 388.9057586548166, 388.64729869880443, 389.63139958072384, 390.1488269998545, 388.91138544610874, 389.57905936367644, 389.4504769524226, 389.22625067001826, 389.79318503855836, 389.6127033215785, 389.1037017285112, 389.73580729121386, 389.31584083722316, 389.40881680859326, 389.323138549745, 389.2451419219243, 389.3922639055792, 388.796221854985, 388.97157341770327, 388.99869472269773, 389.07910256875385, 389.31098859776864, 388.7419881503485, 388.3143923895009, 389.0092636036939, 388.229773876986, 388.3854174264933, 388.05911875263905, 387.2643574001875, 387.3554637363693, 387.39402352923776, 387.02865949856425, 386.9609406993706, 386.7226638865326, 386.7789454376136, 386.7732373979207, 385.4992193760886, 385.5215522567318, 384.82993023752476, 384.24706703575663, 384.2926857786712, 383.6152585260428, 382.8533589897697, 384.67220539781823, 383.42475883749, 382.5620400263505, 376.3078422852331};
	
	int loop1 = s_cont.size();
	int loop2 = 36;

	for(int b=1; b<loop2; b++){
		for(int i=0; i<loop1; i++){
			s_cont.push_back(s_cont[i]);
		}
	}

	double c = 0.1;
	s_cont.insert(s_cont.begin(),c);

	if(s_cont.size() != vectorscin.size()){cout<<"ERROR in calibration!"<<endl;}

	std::vector<double> Calib_vectorScin;

	for(uint i=0; i<vectorscin.size(); i++){
		Calib_vectorScin.push_back(vectorscin[i]*(1.0/s_cont[i]));
	}

	return Calib_vectorScin;
}

std::vector<double> calibcher(std::vector<double> vectorcher){
	std::vector<double> c_cont;
	c_cont = {99.869475119436, 99.5892930698953, 99.36819945465197, 99.32851833492266, 99.23212488346117, 99.21702593539936, 99.03028402205433, 99.15912908522857, 99.2060316999247, 99.30522580785956, 99.11344505994408, 99.21283817184117, 99.25155522412936, 99.29526075724375, 99.32726514547564, 99.21535433745935, 99.24785175107958, 99.2500389422396, 99.2947777840262, 99.20923521702923, 99.30797870907622, 99.31656510767378, 99.21947678681387, 99.32713686590662, 99.32851884709454, 99.305911407397, 99.30508949248147, 99.34913871640343, 99.19261358390726, 99.16893956672706, 99.26948077177317, 99.31331547462389, 99.30086159029355, 99.37873056276479, 99.35080907470335, 99.42847084872328, 99.40510235111934, 99.30670688191663, 99.38448289299362, 99.43075203396795, 99.34555944179634, 99.38032213687576, 99.3653611014606, 99.33974923759105, 99.36938995243544, 99.41133031804752, 99.33192563010871, 99.18049429253469, 99.39705625808192, 99.24849638998961, 99.2286007870197, 99.09928229987484, 99.0993256503711, 99.19364178823666, 99.0555902591915, 99.0621809829577, 98.97675936409205, 99.13401462382535, 99.03520349760886, 99.01410173081639, 98.805439062685, 98.8931348280415, 98.84927421673683, 98.71812808136468, 98.76605765544863, 98.8150716743538, 98.664618195543, 98.52979814367153, 98.39802169218106, 98.48346797819356, 98.46914717338315, 98.21105203135059, 98.38198816987828, 98.04422762750505, 96.54297571859878};

	int loop1 = c_cont.size();
	int loop2 = 36;

	for(int b=1; b<loop2; b++){
		for(int i=0; i<loop1; i++){
			c_cont.push_back(c_cont[i]);
		}
	}

	double c = 0.1;
	c_cont.insert(c_cont.begin(),c);

	if(c_cont.size() != vectorcher.size()){cout<<"ERROR in calibration!"<<endl;}

	std::vector<double> Calib_vectorCher;

	for(uint i=0; i<vectorcher.size(); i++){
		Calib_vectorCher.push_back(vectorcher[i]*(1.0/c_cont[i]));
	}

	return Calib_vectorCher;
}

************************** new calibration constants *********************/



std::vector<double> calibscin(std::vector<double> vectorscin){
	std::vector<double> s_cont;
	s_cont = {408.21638950554075, 408.3954472740771, 407.1870232421094, 406.63875945884087, 404.8060585388971, 403.97304819147996, 403.3691105878475, 403.49367909804056, 404.55647780600043, 405.58591491094637, 403.9575182245898, 404.4757730162475, 404.72249522199195, 405.272159576985, 404.74332809708255, 404.83205898107536, 405.23195412471205, 404.9766105533868, 404.9085068798063, 404.9314555180952, 404.67532710488985, 404.58364980855805, 405.012793566413, 405.0007315500301, 404.30902206187204, 405.6974274788762, 405.2261341502687, 405.63975175649347, 404.90683641527, 404.37034541526305, 405.67260217215875, 405.5109490861691, 404.2898135363692, 405.07073526391474, 405.58981257625425, 405.3751447994642, 405.36549518339785, 405.3332161707569, 404.88956759976287, 405.37027184803094, 404.8980725551248, 405.34774082392767, 405.2984093045488, 405.14372480308344, 405.19187487160525, 405.03757034167137, 405.16280927227615, 404.7829216539207, 405.03107640207867, 404.7292557576276, 404.8025372723253, 403.9177916263665, 404.7460239584375, 403.96821450150077, 404.1905949169899, 404.1704924951662, 403.16496315846314, 402.2360298379118, 403.3863719919289, 402.9762332238292, 403.15699339382735, 403.4020052256797, 402.3032561236677, 402.8453577277423, 401.11356268338346, 401.3504783424065, 400.94087925309395, 400.29569405733, 400.0328154316862, 399.5130445431503, 398.66148407548866, 399.83880015591535, 398.96289406538807, 398.42261837089694, 391.76612693948175};

	
	int loop1 = s_cont.size();
	int loop2 = 36;

	for(int b=1; b<loop2; b++){
		for(int i=0; i<loop1; i++){
			s_cont.push_back(s_cont[i]);
		}
	}

	double c = 0.1;
	s_cont.insert(s_cont.begin(),c);

	if(s_cont.size() != vectorscin.size()){cout<<"ERROR in calibration!"<<endl;}

	std::vector<double> Calib_vectorScin;

	for(uint i=0; i<vectorscin.size(); i++){
		Calib_vectorScin.push_back(vectorscin[i]*(1.0/s_cont[i]));
	}

	return Calib_vectorScin;
}

std::vector<double> calibcher(std::vector<double> vectorcher){
	std::vector<double> c_cont;
	c_cont = {103.08779161895677, 102.91302749597065, 102.69865952763615, 102.61869191270468, 102.54928716539662, 102.48068194031679, 102.49984890080964, 102.35556540203991, 102.47969263317724, 102.6281510005559, 102.43322742473204, 102.47810836409134, 102.55371034296142, 102.67118096060427, 102.67297232291142, 102.48284061965019, 102.5649981010228, 102.56155933915096, 102.67809243921879, 102.56067521092992, 102.60224889784466, 102.63726587197354, 102.63191774143888, 102.76496337880408, 102.6929637252195, 102.60491403169074, 102.85913301772406, 102.741217657914, 102.69546934772463, 102.67035622618218, 102.69304228926421, 102.75886941001674, 102.75976221892324, 102.731492956408, 102.7188845221274, 102.77429845330465, 102.78649420797491, 102.75140309520445, 102.70051794706535, 102.68996042906552, 102.78365100098196, 102.8153738834064, 102.71292597825087, 102.73146416207084, 102.6450394621172, 102.61404003462839, 102.66675609739092, 102.60991640602225, 102.750246685674, 102.62575682868824, 102.42720794074478, 102.51305416968992, 102.52098979376447, 102.59751750679058, 102.45780037787654, 102.53083482963227, 102.47068539942974, 102.5721049950492, 102.56599170316093, 102.46469174495641, 102.19238017547394, 102.28148980648412, 102.19817435184497, 102.1330715125064, 102.09230341456059, 102.05765775486448, 101.9644426420847, 101.96014956820567, 101.85273676485993, 101.93311307596035, 101.96637882465569, 101.68716060542853, 101.55050000833062, 101.67603040894112, 99.77195006099979};

	int loop1 = c_cont.size();
	int loop2 = 36;

	for(int b=1; b<loop2; b++){
		for(int i=0; i<loop1; i++){
			c_cont.push_back(c_cont[i]);
		}
	}

	double c = 0.1;
	c_cont.insert(c_cont.begin(),c);

	if(c_cont.size() != vectorcher.size()){cout<<"ERROR in calibration!"<<endl;}

	std::vector<double> Calib_vectorCher;

	for(uint i=0; i<vectorcher.size(); i++){
		Calib_vectorCher.push_back(vectorcher[i]*(1.0/c_cont[i]));
	}

	return Calib_vectorCher;
}


void recalibrate(std::vector<fastjet::PseudoJet> * jet, TProfile * h_cal)
{

  for (unsigned int i = 0; i < jet->size(); ++i){
    double eshare = energyShare(jet->at(i));
    double calConst = h_cal->GetBinContent(h_cal->FindBin(eshare));
    if (calConst != 0) calConst = 1./calConst;
    else calConst = 1;
    jet->at(i).reset_momentum(calConst * jet->at(i).px(),
			      calConst * jet->at(i).py(),
			      calConst * jet->at(i).pz(),
			      calConst * jet->at(i).E());
			      
  }

}  
/*********************************************************/
