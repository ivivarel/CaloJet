#! /usr/bin/env python 

import ROOT 

ROOT.gStyle.SetOptFit()

ifile = ROOT.TFile.Open("bonsai.bonsai.testfile.root")
t = ifile.Get("MyTree")

c = ROOT.TCanvas()

# make histograms of total response

h_rec = ROOT.TH1F("h_mass","",160,100,180)
t.Draw("higgsrecem_m >> h_mass", "j1t_theta > 0.7 && j1t_theta < 2.4 && j2t_theta > 0.7 && j2t_theta < 2.4 && j1t_E > 10 && j1s_E > 1");
h_rec.SetXTitle("m_{jj} [GeV]")
h_rec.Draw()
h_rec.Fit("gaus","","",120,140)
c.Print("RecoHiggsMass.pdf")

h_scin = ROOT.TH1F("h_mass_scin","",160,100,180)
t.Draw("higgsscint_m >> h_mass_scin","j1t_theta > 0.7 && j1t_theta < 2.4 && j2t_theta > 0.7 && j2t_theta < 2.4 && j1t_E > 10 && j1s_E > 1");
h_cher = ROOT.TH1F("h_mass_cher","",160,100,180)
t.Draw("higgscher_m >> h_mass_cher","j1t_theta > 0.7 && j1t_theta < 2.4 && j2t_theta > 0.7 && j2t_theta < 2.4 && j1t_E > 10 && j1s_E > 1");

h_rec.Draw()
h_scin.SetLineColor(ROOT.kRed)
h_cher.SetLineColor(ROOT.kGreen)
h_rec.Draw()
h_cher.Draw("same")
h_scin.Draw("same")

h_leg = ROOT.TLegend(0.15,0.6,0.35,0.8)
h_leg.SetShadowColor(10)
h_leg.SetLineColor(10)
h_leg.AddEntry(h_rec,"Combined","l")
h_leg.AddEntry(h_scin,"S-channel","l")
h_leg.AddEntry(h_cher,"C-channel","l")
h_leg.Draw()

c.Print("CompareChannels.pdf")


### make histograms to understand the response

h_col_scint_R = ROOT.TH2F("h_col_scint_R", "", 32,0,3.2,40,0,2)
t.Draw("j1s_E/j1t_E:j1t_theta >> h_col_scint_R","j1t_E > 10 && j1s_E > 10","COLZ")
h_col_scint_R.SetXTitle("#theta_{#gamma}")
h_col_scint_R.SetYTitle("E_{cher}/#E_{#gamma}")
h_col_scint_R.Draw("COLZ")
c.Print("ScintResponse2D.pdf")

h_col_cher_R = ROOT.TH2F("h_col_cher_R", "", 32,0,3.2,40,0,2)
t.Draw("j1c_E/j1t_E:j1t_theta >> h_col_cher_R","j1t_E > 10 && j1c_E > 10","COLZ")
h_col_cher_R.SetXTitle("#theta_{#gamma}")
h_col_cher_R.SetYTitle("E_{cher}/#E_{#gamma}")
h_col_cher_R.Draw("COLZ")
c.Print("CherResponse2D.pdf")

h_prof_scint_R = ROOT.TProfile("h_prof_scint_R", "", 64,0,3.2)
t.Draw("j1s_E/j1t_E:j1t_theta >> h_prof_scint_R","j1t_E > 10 && j1s_E > 10","PROF")
h_prof_scint_R.SetXTitle("#theta_{#gamma}")
h_prof_scint_R.SetYTitle("E_{cher}/E_{#gamma}")
h_prof_scint_R.SetMinimum(1)
h_prof_scint_R.SetMaximum(1.1)
h_prof_scint_R.Draw()
c.Print("ScintResponseProf.pdf")

h_prof_cher_R = ROOT.TProfile("h_prof_cher_R", "", 64,0,3.2)
t.Draw("j1c_E/j1t_E:j1t_theta >> h_prof_cher_R","j1t_E > 10 && j1c_E > 10","PROF")
h_prof_cher_R.SetXTitle("#theta_{#gamma}")
h_prof_cher_R.SetYTitle("E_{cher}/E_{#gamma}")
h_prof_cher_R.SetMinimum(1)
h_prof_cher_R.SetMaximum(1.1)
h_prof_cher_R.Draw()
c.Print("CherResponseProf.pdf")

################ Plot hybrid response

h_rec_def = ROOT.TH1F("h_mass_def","",160,100,180)
t.Draw("higgsrecem_m >> h_mass_def", "j1t_theta > 0.7 && j1t_theta < 2.4 && j2t_theta > 0.7 && j2t_theta < 2.4 && j1t_E > 10 && j1s_E > 1");

h_hyb_ene = ROOT.TH1F("h_mass_hyb_ene","",160,100,180)
t.Draw("higgsHybErec_m >> h_mass_hyb_ene", "j1t_theta > 0.7 && j1t_theta < 2.4 && j2t_theta > 0.7 && j2t_theta < 2.4 && j1t_E > 10 && j1s_E > 10");
h_hyb_ene.SetXTitle("m_{jj} [GeV]")

h_hyb_dir = ROOT.TH1F("h_mass_hyb_dir","",160,100,180)
t.Draw("higgsHybDirRec_m >> h_mass_hyb_dir","j1t_theta > 0.7 && j1t_theta < 2.4 && j2t_theta > 0.7 && j2t_theta < 2.4 && j1t_E > 10 && j1s_E > 10");



h_hyb_ene.SetLineColor(ROOT.kCyan)
h_hyb_dir.SetLineColor(ROOT.kGreen)
h_hyb_dir.Fit("gaus","","",123,126)
h_hyb_dir.Draw()
h_rec_def.Draw("same")
h_hyb_ene.Draw("same")

h_leg = ROOT.TLegend(0.15,0.6,0.35,0.8)
h_leg.SetShadowColor(10)
h_leg.SetLineColor(10)
h_leg.AddEntry(h_rec,"Combined","l")
h_leg.AddEntry(h_hyb_dir,"Using truth energy","l")
h_leg.AddEntry(h_hyb_ene,"Using truth direction","l")
h_leg.Draw()

c.Print("HybridPlots.pdf")
