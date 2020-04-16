#! /usr/bin/env python 

import ROOT 

ROOT.gStyle.SetOptFit()

doHybridPlots = True
doHgamgamCal = True
doRecoHiggsMass = True
doResponse = True
doCompareCalib = True

ifile = ROOT.TFile.Open("bonsai.testfile.root")
t = ifile.Get("MyTree")

c = ROOT.TCanvas()

# make histograms of total response

if doRecoHiggsMass:

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

if doResponse: 
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
    h_prof_scint_R.SetMinimum(0.95)
    h_prof_scint_R.SetMaximum(1.05)
    h_prof_scint_R.Draw()
    c.Print("ScintResponseProf.pdf")

    h_prof_cher_R = ROOT.TProfile("h_prof_cher_R", "", 64,0,3.2)
    t.Draw("j1c_E/j1t_E:j1t_theta >> h_prof_cher_R","j1t_E > 10 && j1c_E > 10","PROF")
    h_prof_cher_R.SetXTitle("#theta_{#gamma}")
    h_prof_cher_R.SetYTitle("E_{cher}/E_{#gamma}")
    h_prof_cher_R.SetMinimum(0.95)
    h_prof_cher_R.SetMaximum(1.05)
    h_prof_cher_R.Draw()
    c.Print("CherResponseProf.pdf")

################ Plot hybrid response

if doHybridPlots:

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

############## Compute calibration correction constants

if doHgamgamCal:

    h_calib_scint1 = ROOT.TProfile("h_calib_scint1","",40,0,1)
    t.Draw("j1s_E/j1t_E:j1s_eshare >> h_calib_scint1","j1t_E > 10 && j1t_theta < 2.4 && j1t_theta > 0.7","PROF")
    h_calib_scint1.SetXTitle("E_{tower 1}/E_{tower 2}")
    
    h_calib_scint2 = ROOT.TProfile("h_calib_scint2","",40,0,1)
    t.Draw("j2s_E/j2t_E:j2s_eshare >> h_calib_scint2","j2t_E > 10 && j2t_theta < 2.4 && j2t_theta > 0.7","PROF")
    h_calib_scint2.SetXTitle("E_{tower 1}/E_{tower 2}")
    
    h_calib_cher1 = ROOT.TProfile("h_calib_cher1","",40,0,1)
    t.Draw("j1c_E/j1t_E:j1c_eshare >> h_calib_cher1","j1t_E > 10 && j1t_theta < 2.4 && j1t_theta > 0.7","PROF")
    h_calib_cher1.SetXTitle("E_{tower 1}/E_{tower 2}")
    
    h_calib_cher2 = ROOT.TProfile("h_calib_cher2","",40,0,1)
    t.Draw("j2c_E/j2t_E:j2c_eshare >> h_calib_cher2","j2t_E > 10 && j2t_theta < 2.4 && j2t_theta > 0.7","PROF")
    h_calib_cher2.SetXTitle("E_{tower 1}/E_{tower 2}")
    
    calfile = ROOT.TFile("calfile_DR_hgamgam_check.root","recreate")
    calfile.cd()
    h_calib_scint1.Write()
    h_calib_scint2.Write()
    h_calib_cher1.Write()
    h_calib_cher2.Write()
    calfile.Close()

if doCompareCalib:
    ROOT.gStyle.SetOptStat(0)
    calfile = ROOT.TFile("calfile_DR_hgamgam_check.root")
    if not calfile.IsOpen():
        print "Cannot find calibration file"
    else:
        h_calib_scint1 = calfile.Get("h_calib_scint1")
        h_calib_scint2 = calfile.Get("h_calib_scint2")
        h_calib_cher1 = calfile.Get("h_calib_cher1")
        h_calib_cher2 = calfile.Get("h_calib_cher2")
        h_calib_scint1.SetMarkerStyle(20)
        h_calib_scint2.SetMarkerStyle(23)
        h_calib_cher1.SetMarkerStyle(24)
        h_calib_cher2.SetMarkerStyle(32)
        h_axes = ROOT.TH2F("h_axes","",100,0,1.05,100,0.7,1.2)
        h_axes.SetXTitle("E_{tower 1}/E_{tower 2}")
        h_axes.SetYTitle("E_{jet}/E_{\gamma}")
        h_axes.Draw()
        h_calib_scint1.Draw("same")
        h_calib_scint2.Draw("same")
        h_calib_cher1.Draw("same")
        h_calib_cher2.Draw("same")

        h_leg = ROOT.TLegend(0.6,0.6,0.8,0.8)                                                                                  

        h_leg.SetShadowColor(10)                                                                         
        h_leg.SetLineColor(10)                                                                    
        h_leg.AddEntry(h_calib_scint1,"E_{scin1}","p")                                                                        
        h_leg.AddEntry(h_calib_scint2,"E_{scin2}","p")                                                                         
        h_leg.AddEntry(h_calib_cher1,"E_{cher1}","p")                                                                      
        h_leg.AddEntry(h_calib_cher2,"E_{cher2}","p")                                                                                      
        h_leg.Draw()                     
        
        c.Print("CompareCalib.pdf")
