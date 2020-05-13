#! /usr/bin/env python 

import ROOT,argparse

def CopyTreeToFile(origin_file = None, dest_file = None, origin_treename = "truth"):
    print "Copying tree " + origin_treename + " from file " + origin_file + " to " + dest_file
    # open the files if they are there
    f1 = None
    f2 = None 
    try:
        f1 = ROOT.TFile.Open(origin_file)
    except: 
        print "Cannot find file " + origin_file 
        return False 
    try: 
        f2 = ROOT.TFile.Open(dest_file,"update")
    except: 
        print "Cannot find file " + dest_file 
        return False
        
    t1 = f1.Get(origin_treename)
    
    f2.cd()
    t3 = t1.CloneTree()
    t3.Write()
    f2.Close()
    f1.Close()    

def main():
    print "in main"
    import glob 
#    for myfile in glob.glob("/its/home/iv41/scratch/G4/Collisions/Hgamgam/root/run_0001_Hgamgam_*_truth.root"):
    for myfile in glob.glob("/its/home/iv41/scratch/G4/Collisions/jet_highstat_leakage_B/hgamgam/hgamgam50k_*truth.root"):
        infile = myfile 
        outfile = myfile.split("_truth")[0] + ".root"
        CopyTreeToFile(origin_file = infile, dest_file = outfile)

if __name__ == "__main__":
                            main()         
