   Double_t        Energyem;
   Double_t        EnergyScin;
   Double_t        EnergyCher;
   Double_t        neutrinoleakage;
   Double_t        leakage;
   Double_t        NofCherenkovDetected;
   Double_t        EnergyTot;
   Double_t        PrimaryParticleEnergy;
   Char_t          PrimaryParticleName[13];
   vector<double>  *VectorSignalsR;
   vector<double>  *VectorSignalsL;
   vector<double>  *VectorSignalsCherR;
   vector<double>  *VectorSignalsCherL;
   vector<double>  *VectorL;
   vector<double>  *VectorR;
   vector<double>  *VectorL_loop;
   vector<double>  *VectorR_loop;

   // List of branches
   TBranch        *b_Energyem;   //!
   TBranch        *b_EnergyScin;   //!
   TBranch        *b_EnergyCher;   //!
   TBranch        *b_neutrinoleakage;   //!
   TBranch        *b_leakage;   //!
   TBranch        *b_NofCherenkovDetected;   //!
   TBranch        *b_EnergyTot;   //!
   TBranch        *b_PrimaryParticleEnergy;   //!
   TBranch        *b_PrimaryParticleName;   //!
   TBranch        *b_VectorSignalsR;   //!
   TBranch        *b_VectorSignalsL;   //!
   TBranch        *b_VectorSignalsCherR;   //!
   TBranch        *b_VectorSignalsCherL;   //!
   TBranch        *b_VectorL;   //!
   TBranch        *b_VectorR;   //!
   TBranch        *b_VectorL_loop;
   TBranch        *b_VectorR_loop;
   UInt_t          mcs_n;
   vector<float>   *mcs_E;
   vector<float>   *mcs_pt;
   vector<float>   *mcs_m;
   vector<float>   *mcs_eta;
   vector<float>   *mcs_phi;
   vector<int>     *mcs_status;
   vector<int>     *mcs_barcode;
   vector<int>     *mcs_pdgId;
   vector<float>   *mcs_charge;
   vector<float>   *mcs_vx_x;
   vector<float>   *mcs_vx_y;
   vector<float>   *mcs_vx_z;

   // List of branches
   TBranch        *b_mcs_n;   //!
   TBranch        *b_mcs_E;   //!
   TBranch        *b_mcs_pt;   //!
   TBranch        *b_mcs_m;   //!
   TBranch        *b_mcs_eta;   //!
   TBranch        *b_mcs_phi;   //!
   TBranch        *b_mcs_status;   //!
   TBranch        *b_mcs_barcode;   //!
   TBranch        *b_mcs_pdgId;   //!
   TBranch        *b_mcs_charge;   //!
   TBranch        *b_mcs_vx_x;   //!
   TBranch        *b_mcs_vx_y;   //!
   TBranch        *b_mcs_vx_z;   //!
