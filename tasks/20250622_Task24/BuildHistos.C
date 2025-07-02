/*
Author: Giulia Cossutti

Macro to fill histograms for signal and backgrounds
2 categories of events:
- A: Events with 1 Large Jet
- B: Events with >=2 Large Jets

From inside the /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes directory run with
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250622_Task24/BuildHistos.C'("signal_EJ")' &
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250622_Task24/BuildHistos.C'("bkg_Zjets_q")' &
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250622_Task24/BuildHistos.C'("bkg_Zjets_g")' &
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250622_Task24/BuildHistos.C'("bkg_ZZ")' &
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250622_Task24/BuildHistos.C'("bkg_HZ_SM")' &

Pre-selection on events:
- At least 2 small jets
- 2 charged leptons with same flavour, opposite charge, | Invmass - mZ | < tolerance = 10 GeV, ZpT > 20 GeV

Selection on events:
- >= 1 Large Jet (to fit in category A or B)
- number of b-tagged small jets >= 1 (for both A and B events)
- A: DeltaR between first 2 small jets < 2.2
- B: DeltaR between first 2 small jets < 4.8
- A: pT of sum of 2 leptons by Z decay > 96 GeV/c
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "/gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/classes/DelphesClasses.h"
#include "/gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#endif

//------------------------------------------------------------------------------
//----------- FUNCTIONS -----------------------
//===============================
// Function to compute DeltaR between 2 electrons
Double_t DeltaR(Electron* el1, Electron* el2){
 return sqrt( pow(el1->Eta - el2->Eta,2) + pow(el1->Phi - el2->Phi,2) );
}
//===============================
// Function to compute DeltaR between 2 muons
Double_t DeltaR(Muon* mu1, Muon* mu2){
 return sqrt( pow(mu1->Eta - mu2->Eta,2) + pow(mu1->Phi - mu2->Phi,2) );
}
//===============================
// Function to compute DeltaR between 2 particles
Double_t DeltaR(GenParticle* part1, GenParticle* part2){
 return sqrt( pow(part1->Eta - part2->Eta,2) + pow(part1->Phi - part2->Phi,2) );
}
//===============================
// Function to compute DeltaR between 2 jets
Double_t DeltaR(Jet* jet1, Jet* jet2){
 return sqrt( pow(jet1->Eta - jet2->Eta,2) + pow(jet1->Phi - jet2->Phi,2) );
}
//===============================
// Function to compute DeltaR between a jet and a track
Double_t DeltaR(Jet* jet, Track* trk){
 return sqrt( pow(jet->Eta - trk->Eta,2) + pow(jet->Phi - trk->Phi,2) );
}
//===============================
// Function to compute DeltaR between 2 tracks
Double_t DeltaR(Track* trk1, Track* trk2){
 return sqrt( pow(trk1->Eta - trk2->Eta,2) + pow(trk1->Phi - trk2->Phi,2) );
}
//===============================
// Function to compute DeltaEta between 2 jets
Double_t DeltaEta(Jet* jet1, Jet* jet2){
 return jet1->Eta - jet2->Eta;
}
//===============================
// Function to compute DeltaPhi between 2 jets
Double_t DeltaPhi(Jet* jet1, Jet* jet2){
 return jet1->Phi - jet2->Phi;
}
//===============================
// Function to calculate Prompt Track Fraction without z cut
Double_t PTFnoz(ExRootTreeReader* treeReader,Int_t entry,TClonesArray *branchTrackPTF, Jet* jet){
  TClonesArray *branchTrack = branchTrackPTF;

  // Load selected branches with data from specified event
  treeReader->ReadEntry(entry);

  // Sum of PT of selected tracks
  Double_t sumPT = 0.0;

  // Number of used tracks
  Int_t count = 0;

  // Loop on tracks
  for (Int_t it = 0; it < branchTrack->GetEntries(); it++){
    Track* trk = (Track*)branchTrack->At(it);

    // Track definition for Large Radius Tracking algorithm
    if( (trk->D0 >= 300) || (trk->DZ >= 500) || (trk->PT <= 1) || (abs(trk->Eta) >= 3.0) ) continue;

    // Track within DeltaR = 1.2 from the jet
    if( DeltaR(jet, trk) >= 1.2 ) continue;

    // Cut on d0
    if( trk->D0 >= 2.5 * (0.0463/(trk->PT) + 0.0195) ) continue;

    // Increment track counter
    count += 1;

    // Add track PT to sum
    sumPT += trk->PT;
  }

  //cout << "Tracks used in PTFnoz for event " << entry << " : " << count << endl; 

  return sumPT/(jet->PT);
}
//===============================
// Function to calculate Prompt Track Fraction with z cut
Double_t PTF(ExRootTreeReader* treeReader,Int_t entry,TClonesArray *branchTrackPTF, Jet* jet, Double_t zPV){
  TClonesArray *branchTrack = branchTrackPTF;

  // Load selected branches with data from specified event
  treeReader->ReadEntry(entry);

  // Sum of PT of selected tracks
  Double_t sumPT = 0.0;

  // Number of used tracks
  Int_t count = 0;

  // Loop on tracks
  for (Int_t it = 0; it < branchTrack->GetEntries(); it++){
    Track* trk = (Track*)branchTrack->At(it);

    // Track definition for Large Radius Tracking algorithm
    if( (trk->D0 >= 300) || (trk->DZ >= 500) || (trk->PT <= 1) || (abs(trk->Eta) >= 3.0) ) continue;

    // Track within DeltaR = 1.2 from the jet
    if( DeltaR(jet, trk) >= 1.2 ) continue;

    // Cut on d0
    if( trk->D0 >= 2.5 * (0.0463/(trk->PT) + 0.0195) ) continue;

    // Cut on z0
    if( abs(zPV - trk->DZ) >= 10 ) continue;

    // Increment track counter
    count += 1;

    // Add track PT to sum
    sumPT += trk->PT;
  }

  //cout << "Tracks used in PTF for event " << entry << " : " << count << endl; 

  return sumPT/(jet->PT);
}
//===============================
// Function to calculate Prompt Track Fraction with z cut for small jets (DeltaR=0.5)
Double_t PTFsmall(ExRootTreeReader* treeReader,Int_t entry,TClonesArray *branchTrackPTF, Jet* jet, Double_t zPV){
  TClonesArray *branchTrack = branchTrackPTF;

  // Load selected branches with data from specified event
  treeReader->ReadEntry(entry);

  // Sum of PT of selected tracks
  Double_t sumPT = 0.0;

  // Number of used tracks
  Int_t count = 0;

  // Loop on tracks
  for (Int_t it = 0; it < branchTrack->GetEntries(); it++){
    Track* trk = (Track*)branchTrack->At(it);

    // Track definition for Large Radius Tracking algorithm
    if( (trk->D0 >= 300) || (trk->DZ >= 500) || (trk->PT <= 1) || (abs(trk->Eta) >= 3.0) ) continue;

    // Track within DeltaR = 0.5 from the jet
    if( DeltaR(jet, trk) >= 0.5 ) continue;

    // Cut on d0
    if( trk->D0 >= 2.5 * (0.0463/(trk->PT) + 0.0195) ) continue;

    // Cut on z0
    if( abs(zPV - trk->DZ) >= 10 ) continue;

    // Increment track counter
    count += 1;

    // Add track PT to sum
    sumPT += trk->PT;
  }

  //cout << "Tracks used in PTF for event " << entry << " : " << count << endl; 

  return sumPT/(jet->PT);
}
//===============================
// Function to calculate 2-point Energy Correlation Function of a large jet
Double_t ECF2(ExRootTreeReader* treeReader,Int_t entry,TClonesArray *branchTrackECF, Jet* jet){
  TClonesArray *branchTrack = branchTrackECF;

  // Load selected branches with data from specified event
  treeReader->ReadEntry(entry);

  // ECF2
  Double_t ECF2 = 0.0;

  // Loop1 on tracks
  for (Int_t i = 0; i < branchTrack->GetEntries(); i++){
    Track* trk1 = (Track*)branchTrack->At(i);

    // Track definition for Large Radius Tracking algorithm
    if( (trk1->D0 >= 300) || (trk1->DZ >= 500) || (trk1->PT <= 1) || (abs(trk1->Eta) >= 3.0) ) continue;

    // Track within DeltaR = 1.0 from the jet
    if( DeltaR(jet, trk1) >= 1.0 ) continue;

    // Loop2 on tracks
    for (Int_t j = 0; j < i; j++){
      Track* trk2 = (Track*)branchTrack->At(j);

      // Track definition for Large Radius Tracking algorithm
      if( (trk2->D0 >= 300) || (trk2->DZ >= 500) || (trk2->PT <= 1) || (abs(trk2->Eta) >= 3.0) ) continue;

      // Track within DeltaR = 1.0 from the jet
      if( DeltaR(jet, trk2) >= 1.0 ) continue;

      // Add to ECF2
      ECF2 += trk1->PT * trk2->PT * DeltaR(trk1, trk2);
    }
  }

  return ECF2;
}
//===============================
// Function to calculate 2-point Energy Correlation Function of a small jet
Double_t ECF2small(ExRootTreeReader* treeReader,Int_t entry,TClonesArray *branchTrackECF, Jet* jet){
  TClonesArray *branchTrack = branchTrackECF;

  // Load selected branches with data from specified event
  treeReader->ReadEntry(entry);

  // ECF2
  Double_t ECF2 = 0.0;

  // Loop1 on tracks
  for (Int_t i = 0; i < branchTrack->GetEntries(); i++){
    Track* trk1 = (Track*)branchTrack->At(i);

    // Track definition for Large Radius Tracking algorithm
    if( (trk1->D0 >= 300) || (trk1->DZ >= 500) || (trk1->PT <= 1) || (abs(trk1->Eta) >= 3.0) ) continue;

    // Track within DeltaR = 0.4 from the jet
    if( DeltaR(jet, trk1) >= 0.4 ) continue;

    // Loop2 on tracks
    for (Int_t j = 0; j < i; j++){
      Track* trk2 = (Track*)branchTrack->At(j);

      // Track definition for Large Radius Tracking algorithm
      if( (trk2->D0 >= 300) || (trk2->DZ >= 500) || (trk2->PT <= 1) || (abs(trk2->Eta) >= 3.0) ) continue;

      // Track within DeltaR = 0.4 from the jet
      if( DeltaR(jet, trk2) >= 0.4 ) continue;

      // Add to ECF2
      ECF2 += trk1->PT * trk2->PT * DeltaR(trk1, trk2);
    }
  }

  return ECF2;
}
//===============================
//----------- END FUNCTIONS ---------------------


void BuildHistos(const char *inFile)
{

  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");

  TString inname;
  for (Int_t i = 1; i <= 10; i++){
    inname.Form("/gfsvol01/atlas/giuliac/condor/20250521_Task20/root_outputs/50K_%s_%d.root",inFile,i);
    chain.Add(inname);
    cout << "File " << inname << " added to chain." << endl;
  }

  if( (std::strcmp(inFile,"bkg_Zjets_q") == 0) || (std::strcmp(inFile,"bkg_Zjets_g") == 0) ){
    for (Int_t i = 11; i <= 50; i++){
      inname.Form("/eos/infnts/atlas/giuliac/condor/20250527_Task21/root_outputs/50K_%s_%d.root",inFile,i);
      chain.Add(inname);
      cout << "File " << inname << " added to chain." << endl;
    }
  }

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchReclusteredJet = treeReader->UseBranch("ReclusteredJet");
  TClonesArray *branchReclusteredGenJet = treeReader->UseBranch("ReclusteredGenJet");
  
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
    
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchGenMissingET = treeReader->UseBranch("GenMissingET");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchScalarHT = treeReader->UseBranch("ScalarHT");
  
  TClonesArray *branchJetEnergyScaledJet = treeReader->UseBranch("JetEnergyScaledJet");
  TClonesArray *branchECalTower = treeReader->UseBranch("ECalTower");
  TClonesArray *branchHCalTower = treeReader->UseBranch("HCalTower");
  TClonesArray *branchMuonMomentumSmearing = treeReader->UseBranch("MuonMomentumSmearing");
  
  // Event categories
  std::vector<TString> cat = {"A_","B_"};

  // Book histograms
  // 1D
  std::vector<TH1F*> Invmass_2small;
  std::vector<TH1F*> Invmass_2large;
  std::vector<TH1F*> Invmass_bsmall;

  std::vector<TH1F*> DeltaR_small;
  std::vector<TH1F*> DeltaPhi_small;
  std::vector<TH1F*> DeltaR_large;
  std::vector<TH1F*> DeltaPhi_large;

  std::vector<TH1F*> JetPT1_small;
  std::vector<TH1F*> JetPT2_small;
  std::vector<TH1F*> JetPT1_large;
  std::vector<TH1F*> JetPT2_large;

  std::vector<TH1F*> NConst1_large;
  std::vector<TH1F*> NConst2_large;
  std::vector<TH1F*> N_small;
  std::vector<TH1F*> N_bsmall;
  std::vector<TH1F*> N_large;

  std::vector<TH1F*> PTF1_large;
  std::vector<TH1F*> PTF2_large;
  std::vector<TH1F*> PTF1_small;
  std::vector<TH1F*> PTF2_small;
  std::vector<TH1F*> PTF1_bsmall;

  std::vector<TH1F*> ECF21_large;
  std::vector<TH1F*> ECF22_large;
  std::vector<TH1F*> ECF21_small;
  std::vector<TH1F*> ECF22_small;

  std::vector<TH1F*> Lep2PT;
  std::vector<TH1F*> LepDeltaR;

  // 2D
  std::vector<TH2F*> IMsmall_DeltaRs;
  std::vector<TH2F*> IMlarge_DeltaRl;

  std::vector<TH2F*> PTF1s_PTF2s;
  std::vector<TH2F*> PTF1s_PTF1l;
  std::vector<TH2F*> PTF1bs_PTF1l;

  for(Int_t i = 0; i < cat.size(); ++i){
    Invmass_2small.push_back(new TH1F(cat.at(i) + "Invariant_mass_of_leading_and_subleading_small_jets", "", 200, 0.0, 800.0));
    Invmass_2large.push_back(new TH1F(cat.at(i) + "Invariant_mass_of_leading_and_subleading_large_jets", "", 200, 0.0, 800.0));
    Invmass_bsmall.push_back(new TH1F(cat.at(i) + "Invariant_mass_of_all_btagged_small_jets", "", 200, 0.0, 800.0));

    DeltaR_small.push_back(new TH1F(cat.at(i) + "DeltaR_between_leading_and_subleading_small_jets", "", 100, 0., 10.));
    DeltaPhi_small.push_back(new TH1F(cat.at(i) + "DeltaPhi_between_leading_and_subleading_small_jets", "", 110, -1., 10.));
    DeltaR_large.push_back(new TH1F(cat.at(i) + "DeltaR_between_leading_and_subleading_large_jets", "", 100, 0., 10.));
    DeltaPhi_large.push_back(new TH1F(cat.at(i) + "DeltaPhi_between_leading_and_subleading_large_jets", "", 110, -1., 10.));

    JetPT1_small.push_back(new TH1F(cat.at(i) + "leading_small_jet_pt", "", 95, 20.0, 400.0));
    JetPT2_small.push_back(new TH1F(cat.at(i) + "subleading_small_jet_pt", "", 95, 20.0, 400.0));
    JetPT1_large.push_back(new TH1F(cat.at(i) + "leading_large_jet_pt", "", 90, 40.0, 400.0));
    JetPT2_large.push_back(new TH1F(cat.at(i) + "subleading_large_jet_pt", "", 90, 40.0, 400.0));
 
    NConst1_large.push_back(new TH1F(cat.at(i) + "number_of_constituents_of_leading_large_jet", "", 17, -1.5, 15.5));
    NConst2_large.push_back(new TH1F(cat.at(i) + "number_of_constituents_of_subleading_large_jet", "", 17, -1.5, 15.5));
    N_small.push_back(new TH1F(cat.at(i) + "number_of_small_jets", "", 17, -1.5, 15.5));
    N_bsmall.push_back(new TH1F(cat.at(i) + "number_of_btagged_small_jets", "", 17, -1.5, 15.5));
    N_large.push_back(new TH1F(cat.at(i) + "number_of_large_jets", "", 17, -1.5, 15.5));

    PTF1_large.push_back(new TH1F(cat.at(i) + "Prompt_Track_Fraction_of_leading_large_jet", "", 50, 0., 1.));
    PTF2_large.push_back(new TH1F(cat.at(i) + "Prompt_Track_Fraction_of_subleading_large_jet", "", 50, 0., 1.));
    PTF1_small.push_back(new TH1F(cat.at(i) + "Prompt_Track_Fraction_of_leading_small_jet", "", 50, 0., 1.));
    PTF2_small.push_back(new TH1F(cat.at(i) + "Prompt_Track_Fraction_of_subleading_small_jet", "", 50, 0., 1.));
    PTF1_bsmall.push_back(new TH1F(cat.at(i) + "Prompt_Track_Fraction_of_leading_btagged_small_jet", "", 50, 0., 1.));

    ECF21_large.push_back(new TH1F(cat.at(i) + "2point_Energy_Correlation_div_pT_of_leading_large_jet", "", 105, -10., 200.));
    ECF22_large.push_back(new TH1F(cat.at(i) + "2point_Energy_Correlation_div_pT_of_subleading_large_jet", "", 105, -10., 200.));
    ECF21_small.push_back(new TH1F(cat.at(i) + "2point_Energy_Correlation_div_pT_of_leading_small_jet", "", 105, -10., 200.));
    ECF22_small.push_back(new TH1F(cat.at(i) + "2point_Energy_Correlation_div_pT_of_subleading_small_jet", "", 105, -10., 200.));

    Lep2PT.push_back(new TH1F(cat.at(i) + "pt_of_first_two_leptons", "", 150, 0.0, 600.));
    LepDeltaR.push_back(new TH1F(cat.at(i) + "DeltaR_between_first_two_leptons", "", 100, -0.5, 9.5));

    IMsmall_DeltaRs.push_back(new TH2F(cat.at(i) + "IMsmall_vs_DeltaRsmall", "", 200, 0.0, 800.0, 100, 0.0, 10.0));
    IMlarge_DeltaRl.push_back(new TH2F(cat.at(i) + "IMlarge_vs_DeltaRlarge", "", 200, 0.0, 800.0, 100, 0.0, 10.0));

    PTF1s_PTF2s.push_back(new TH2F(cat.at(i) + "PTF1small_vs_PTF2small", "", 50, 0.0, 1.0, 50, 0.0, 1.0));
    PTF1s_PTF1l.push_back(new TH2F(cat.at(i) + "PTF1small_vs_PTF1large", "", 50, 0.0, 1.0, 50, 0.0, 1.0));
    PTF1bs_PTF1l.push_back(new TH2F(cat.at(i) + "PTF1bsmall_vs_PTF1large", "", 50, 0.0, 1.0, 50, 0.0, 1.0));
  }

  // Counter of used events
  Int_t NofEv = 0;
  Int_t NofEvA = 0;
  Int_t NofEvB = 0;

  // Tolerance on Invariant mass near Z mass [GeV]
  Double_t tolerance = 10.;

  // z of Hard Scattering Primary Vertex
  Double_t zPV = 0.;

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

     //-------------------- Preselection -------------------------

     // Preselection on number of small jets

     // Counter of small and b-tagged small
     Int_t n_smalljets = 0;
     Int_t n_bsmalljets = 0;
     Int_t index_bsmall1 = 0;
     Int_t index_small1 = 0;
     Int_t index_small2 = 0;
     
     // Loop over all small jets
     for(Int_t jetentry = 0; jetentry < branchJet->GetEntries(); ++jetentry)
     {
       Jet *jet = (Jet*) branchJet->At(jetentry);
       //if( abs(jet->Eta) >= 2.5 ) continue;
       n_smalljets += 1;
       if( jet->BTag == 1 ) n_bsmalljets += 1;

       // Indexes of leading and subleading small jets
       if(n_smalljets == 1){
         index_small1 = jetentry;
       }else if(n_smalljets == 2){
         index_small2 = jetentry;
       }
       if(n_bsmalljets == 1){
         index_bsmall1 = jetentry;
       }
     }//end loop over all small jets

     if( (n_smalljets < 2) ) continue;

     //------------------ Selection: CUT1 number of b-tagged Small Jets --------------------

     if( (n_bsmalljets < 1) ) continue;

     //------------------ End of Selection: CUT1 number of b-tagged Small Jets --------------------

     //------------------ Selection: number of Large Jets --------------------

     // Count number of large jets to determine category of the event

     // Counter of large
     Int_t n_largejets = 0;
     Int_t index_large1 = 0;
     Int_t index_large2 = 0;
     
     // Loop over all large jets
     for(Int_t jetentry = 0; jetentry < branchReclusteredJet->GetEntries(); ++jetentry)
     {
       Jet *jet = (Jet*) branchReclusteredJet->At(jetentry);
       //if( abs(jet->Eta) >= 1.8 ) continue;
       n_largejets += 1;

       // Indexes of leading and subleading large jets
       if(n_largejets == 1){
         index_large1 = jetentry;
       }else if(n_largejets == 2){
         index_large2 = jetentry;
       }
     }//end loop over all large jets

     if( (n_largejets < 1) ) continue;

     //------------------ End of Selection: number of Large Jets --------------------

     //------------------ Selection: CUT2 DeltaR between first 2 Small Jets --------------------

       Jet *smalljet1 = (Jet*) branchJet->At(index_small1);
       Jet *smalljet2 = (Jet*) branchJet->At(index_small2);

       if( (n_largejets == 1) & (DeltaR(smalljet1, smalljet2) >= 2.2) ) continue;
       if( (n_largejets >= 2) & (DeltaR(smalljet1, smalljet2) >= 4.8) ) continue;

     //------------------ End of Selection: CUT2 DeltaR between first 2 Small Jets --------------------


     // Preselection on charged leptons from Z decay

     // Is there acceptable Z -> l+ l- ?
     Bool_t isZ2l = kFALSE; 

     // Is Z -> 2 electrons?
     Bool_t isZel = kFALSE;
 
     // Is Z -> 2 muons?
     Bool_t isZmu = kFALSE; 

     // Counter of electrons and muons
     Int_t n_e = 0;
     Int_t n_mu = 0;

     // Loop over electrons
     for(Int_t elentry = 0; elentry < branchElectron->GetEntries(); ++elentry)
     {
       Electron *el = (Electron*) branchElectron->At(elentry);
       n_e += 1;
       if(elentry == 1) break;
     }//end loop over electrons

     if(n_e == 2){
       // Candidate electrons from Z decay
       Electron *el1 = (Electron*) branchElectron->At(0);
       Electron *el2 = (Electron*) branchElectron->At(1);

       // Different charges
       if(el1->Charge != el2->Charge){

         // 4-momentum of electrons
         TLorentzVector *el1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
         el1P4->SetPtEtaPhiM(el1->PT,el1->Eta,el1->Phi,0.511e-3);
         TLorentzVector *el2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
         el2P4->SetPtEtaPhiM(el2->PT,el2->Eta,el2->Phi,0.511e-3);

         // Selection on pT of Z
         if( (*el1P4 + *el2P4).Pt() > 20 ){
         
           // Invariant mass near Z mass
           if( abs( (*el1P4 + *el2P4).M() - 91.188 ) < tolerance ){ 

             isZ2l = kTRUE;
             isZel = kTRUE;

           }
         }
       }
     }

     if(isZ2l == kFALSE){

       // Loop over muons
       for(Int_t muentry = 0; muentry < branchMuon->GetEntries(); ++muentry)
       {
         Muon *mu = (Muon*) branchMuon->At(muentry);
         n_mu += 1;
         if(muentry == 1) break;
       }//end loop over muons

       if(n_mu == 2){
         // Candidate electrons from Z decay
         Muon *mu1 = (Muon*) branchMuon->At(0);
         Muon *mu2 = (Muon*) branchMuon->At(1);

         // Different charges
         if(mu1->Charge != mu2->Charge){

           // 4-momentum of muons
           TLorentzVector *mu1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
           mu1P4->SetPtEtaPhiM(mu1->PT,mu1->Eta,mu1->Phi,105.66e-3);
           TLorentzVector *mu2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
           mu2P4->SetPtEtaPhiM(mu2->PT,mu2->Eta,mu2->Phi,105.66e-3);

           // Selection on pT of Z
           if( (*mu1P4 + *mu2P4).Pt() > 20 ){ 
       
             // Invariant mass near Z mass
             if( abs( (*mu1P4 + *mu2P4).M() - 91.188 ) < tolerance ){

               isZ2l = kTRUE;
               isZmu = kTRUE;

             }
           }
         }
       }
     }

     if( isZ2l == kFALSE ) continue;

     //-------------------- End of Preselection -------------------------

     //------------------ Selection: CUT3 pT of sum of 2 leptons from Z decay --------------------

     if(n_largejets == 1){

       if( isZel == kTRUE ){
         Electron *el1 = (Electron*) branchElectron->At(0);
         Electron *el2 = (Electron*) branchElectron->At(1);

         // 4-momentum of electrons
         TLorentzVector *el1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
         el1P4->SetPtEtaPhiM(el1->PT,el1->Eta,el1->Phi,0.511e-3);
         TLorentzVector *el2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
         el2P4->SetPtEtaPhiM(el2->PT,el2->Eta,el2->Phi,0.511e-3);

         if( (*el1P4 + *el2P4).Pt() <= 96 ) continue;
       }

       if( isZmu == kTRUE ){
         Muon *mu1 = (Muon*) branchMuon->At(0);
         Muon *mu2 = (Muon*) branchMuon->At(1);

         // 4-momentum of muons
         TLorentzVector *mu1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
         mu1P4->SetPtEtaPhiM(mu1->PT,mu1->Eta,mu1->Phi,105.66e-3);
         TLorentzVector *mu2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
         mu2P4->SetPtEtaPhiM(mu2->PT,mu2->Eta,mu2->Phi,105.66e-3);

         if( (*mu1P4 + *mu2P4).Pt() <= 96 ) continue;
       }

     }
     //------------------ End of Selection: CUT3 pT of sum of 2 leptons from Z decay --------------------


     // Increment counter of events passing the selection
     NofEv += 1;
     
     if(n_largejets == 1){
       NofEvA += 1;
     }else if(n_largejets >= 2){
       NofEvB += 1;
     }

     //-------------- End of Pre-selection and Selection -------------------

     // Fill histograms

     // Histograms about leptons
     if( isZel == kTRUE ){
       Electron *el1 = (Electron*) branchElectron->At(0);
       Electron *el2 = (Electron*) branchElectron->At(1);

       // 4-momentum of electrons
       TLorentzVector *el1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       el1P4->SetPtEtaPhiM(el1->PT,el1->Eta,el1->Phi,0.511e-3);
       TLorentzVector *el2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       el2P4->SetPtEtaPhiM(el2->PT,el2->Eta,el2->Phi,0.511e-3);

       if(n_largejets == 1){
         Lep2PT.at(0)->Fill( (*el1P4 + *el2P4).Pt() );
         LepDeltaR.at(0)->Fill( DeltaR(el1, el2) );
       }else if(n_largejets >= 2){
         Lep2PT.at(1)->Fill( (*el1P4 + *el2P4).Pt() );
         LepDeltaR.at(1)->Fill( DeltaR(el1, el2) );
       }
     }

     if( isZmu == kTRUE ){
       Muon *mu1 = (Muon*) branchMuon->At(0);
       Muon *mu2 = (Muon*) branchMuon->At(1);

       // 4-momentum of muons
       TLorentzVector *mu1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       mu1P4->SetPtEtaPhiM(mu1->PT,mu1->Eta,mu1->Phi,105.66e-3);
       TLorentzVector *mu2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       mu2P4->SetPtEtaPhiM(mu2->PT,mu2->Eta,mu2->Phi,105.66e-3);

       if(n_largejets == 1){
         Lep2PT.at(0)->Fill( (*mu1P4 + *mu2P4).Pt() );
         LepDeltaR.at(0)->Fill( DeltaR(mu1, mu2) );
       }else if(n_largejets >= 2){
         Lep2PT.at(1)->Fill( (*mu1P4 + *mu2P4).Pt() );
         LepDeltaR.at(1)->Fill( DeltaR(mu1, mu2) );
       }
     }

     // Histograms about jets

     // Jet definition
     Jet *smallbjet1 = (Jet*) branchJet->At(index_bsmall1);
     Jet *largejet1 = (Jet*) branchReclusteredJet->At(index_large1);

     TLorentzVector *smalljet1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
     smalljet1P4->SetPtEtaPhiM(smalljet1->PT,smalljet1->Eta,smalljet1->Phi,smalljet1->Mass);
     TLorentzVector *smalljet2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
     smalljet2P4->SetPtEtaPhiM(smalljet2->PT,smalljet2->Eta,smalljet2->Phi,smalljet2->Mass);

     TLorentzVector *largejet1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
     largejet1P4->SetPtEtaPhiM(largejet1->PT,largejet1->Eta,largejet1->Phi,largejet1->Mass);


     if(n_bsmalljets >= 2){

       // Calculate InvMass of all btagged small jets
       TLorentzVector alljetbP4 = TLorentzVector(0.0,0.0,0.0,0.0);

       // Loop over small jets
       for(Int_t jetentry = index_bsmall1; jetentry < branchJet->GetEntries(); ++jetentry)
       {
         Jet *jet = (Jet*) branchJet->At(jetentry);
         //if( abs(jet->Eta) >= 2.5 ) continue;
         if( jet->BTag != 1 ) continue;

         TLorentzVector *smalljetbP4 = new TLorentzVector(1.0,1.0,1.0,1.0);
         smalljetbP4->SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
         alljetbP4 += *smalljetbP4;
       }

       if(n_largejets == 1){
         Invmass_bsmall.at(0)->Fill( (alljetbP4).M() );
       }else if(n_largejets >= 2){
         Invmass_bsmall.at(1)->Fill( (alljetbP4).M() );
       }
     }

     if(n_largejets == 1){
       N_small.at(0)->Fill(n_smalljets);
       N_bsmall.at(0)->Fill(n_bsmalljets);
       N_large.at(0)->Fill(n_largejets);
       NConst1_large.at(0)->Fill(largejet1->Constituents.GetEntriesFast());       

       Invmass_2small.at(0)->Fill( (*smalljet1P4 + *smalljet2P4).M()  );

       DeltaR_small.at(0)->Fill( DeltaR(smalljet1, smalljet2) );
       DeltaPhi_small.at(0)->Fill( abs(DeltaPhi(smalljet1, smalljet2)) );

       JetPT1_small.at(0)->Fill( smalljet1->PT );
       JetPT2_small.at(0)->Fill( smalljet2->PT );
       JetPT1_large.at(0)->Fill( largejet1->PT );

       PTF1_small.at(0)->Fill( PTFsmall(treeReader, entry, branchTrack, smalljet1, zPV) );
       PTF2_small.at(0)->Fill( PTFsmall(treeReader, entry, branchTrack, smalljet2, zPV) );
       PTF1_bsmall.at(0)->Fill( PTFsmall(treeReader, entry, branchTrack, smallbjet1, zPV) );
       PTF1_large.at(0)->Fill( PTF(treeReader, entry, branchTrack, largejet1, zPV) );

       ECF21_small.at(0)->Fill( ECF2small(treeReader, entry, branchTrack, smalljet1) / smalljet1->PT );
       ECF22_small.at(0)->Fill( ECF2small(treeReader, entry, branchTrack, smalljet2) / smalljet2->PT );
       ECF21_large.at(0)->Fill( ECF2(treeReader, entry, branchTrack, largejet1) / largejet1->PT );

       IMsmall_DeltaRs.at(0)->Fill( (*smalljet1P4 + *smalljet2P4).M(), DeltaR(smalljet1, smalljet2) );

       PTF1s_PTF2s.at(0)->Fill( PTFsmall(treeReader, entry, branchTrack, smalljet1, zPV) , PTFsmall(treeReader, entry, branchTrack, smalljet2, zPV) );
       PTF1s_PTF1l.at(0)->Fill( PTFsmall(treeReader, entry, branchTrack, smalljet1, zPV) , PTF(treeReader, entry, branchTrack, largejet1, zPV) );
       PTF1bs_PTF1l.at(0)->Fill( PTFsmall(treeReader, entry, branchTrack, smallbjet1, zPV) , PTF(treeReader, entry, branchTrack, largejet1, zPV) );

     }else if(n_largejets >= 2){
       N_small.at(1)->Fill(n_smalljets);
       N_bsmall.at(1)->Fill(n_bsmalljets);
       N_large.at(1)->Fill(n_largejets);
       NConst1_large.at(1)->Fill(largejet1->Constituents.GetEntriesFast());       

       Invmass_2small.at(1)->Fill( (*smalljet1P4 + *smalljet2P4).M()  );

       DeltaR_small.at(1)->Fill( DeltaR(smalljet1, smalljet2) );
       DeltaPhi_small.at(1)->Fill( abs(DeltaPhi(smalljet1, smalljet2)) );

       JetPT1_small.at(1)->Fill( smalljet1->PT );
       JetPT2_small.at(1)->Fill( smalljet2->PT );
       JetPT1_large.at(1)->Fill( largejet1->PT );

       PTF1_small.at(1)->Fill( PTFsmall(treeReader, entry, branchTrack, smalljet1, zPV) );
       PTF2_small.at(1)->Fill( PTFsmall(treeReader, entry, branchTrack, smalljet2, zPV) );
       PTF1_bsmall.at(1)->Fill( PTFsmall(treeReader, entry, branchTrack, smallbjet1, zPV) );
       PTF1_large.at(1)->Fill( PTF(treeReader, entry, branchTrack, largejet1, zPV) );

       ECF21_small.at(1)->Fill( ECF2small(treeReader, entry, branchTrack, smalljet1) / smalljet1->PT );
       ECF22_small.at(1)->Fill( ECF2small(treeReader, entry, branchTrack, smalljet2) / smalljet2->PT );
       ECF21_large.at(1)->Fill( ECF2(treeReader, entry, branchTrack, largejet1) / largejet1->PT );

       IMsmall_DeltaRs.at(1)->Fill( (*smalljet1P4 + *smalljet2P4).M(), DeltaR(smalljet1, smalljet2) );

       PTF1s_PTF2s.at(1)->Fill( PTFsmall(treeReader, entry, branchTrack, smalljet1, zPV) , PTFsmall(treeReader, entry, branchTrack, smalljet2, zPV) );
       PTF1s_PTF1l.at(1)->Fill( PTFsmall(treeReader, entry, branchTrack, smalljet1, zPV) , PTF(treeReader, entry, branchTrack, largejet1, zPV) );
       PTF1bs_PTF1l.at(1)->Fill( PTFsmall(treeReader, entry, branchTrack, smallbjet1, zPV) , PTF(treeReader, entry, branchTrack, largejet1, zPV) );

       // B only histos
       Jet *largejet2 = (Jet*) branchReclusteredJet->At(index_large2);
       TLorentzVector *largejet2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       largejet2P4->SetPtEtaPhiM(largejet2->PT,largejet2->Eta,largejet2->Phi,largejet2->Mass);

       Invmass_2large.at(1)->Fill( (*largejet1P4 + *largejet2P4).M()  );

       DeltaR_large.at(1)->Fill( DeltaR(largejet1, largejet2) );
       DeltaPhi_large.at(1)->Fill( abs(DeltaPhi(largejet1, largejet2)) );

       JetPT2_large.at(1)->Fill( largejet2->PT );

       NConst2_large.at(1)->Fill(largejet2->Constituents.GetEntriesFast());       

       PTF2_large.at(1)->Fill( PTF(treeReader, entry, branchTrack, largejet2, zPV) );

       ECF22_large.at(1)->Fill( ECF2(treeReader, entry, branchTrack, largejet2) / largejet2->PT );

       IMlarge_DeltaRl.at(1)->Fill( (*largejet1P4 + *largejet2P4).M(), DeltaR(largejet1, largejet2) );
     }

  }// end loop over all events

//------------------------------------------------

/*
  cout << "Number of used Events for " << inFile << " :" << endl;
  cout << "Number of Events passing selection: " << NofEv << " out of " << numberOfEntries << " : " << 100.*NofEv/numberOfEntries << " %" << endl;
  cout << "Number of Events of category A passing selection: " << NofEvA << " out of " << numberOfEntries << " : " << 100.*NofEvA/numberOfEntries << " %" << endl;
  cout << "Number of Events of category B passing selection: " << NofEvB << " out of " << numberOfEntries << " : " << 100.*NofEvB/numberOfEntries << " %" << endl;
*/

//-------------------------------------------------------
  // write histograms on file
  TString outname;
  outname.Form("/gfsvol01/atlas/giuliac/plots_and_outputs/20250622_Task24/%s.root",inFile);

  TFile *f = new TFile(outname,"RECREATE");

  for(Int_t i = 0; i < cat.size(); ++i){
    Invmass_2small.at(i)->Write("",TObject::kOverwrite);
    Invmass_2large.at(i)->Write("",TObject::kOverwrite);
    Invmass_bsmall.at(i)->Write("",TObject::kOverwrite);

    DeltaR_small.at(i)->Write("",TObject::kOverwrite);
    DeltaPhi_small.at(i)->Write("",TObject::kOverwrite);
    DeltaR_large.at(i)->Write("",TObject::kOverwrite);
    DeltaPhi_large.at(i)->Write("",TObject::kOverwrite);

    JetPT1_small.at(i)->Write("",TObject::kOverwrite);
    JetPT2_small.at(i)->Write("",TObject::kOverwrite);
    JetPT1_large.at(i)->Write("",TObject::kOverwrite);
    JetPT2_large.at(i)->Write("",TObject::kOverwrite);

    NConst1_large.at(i)->Write("",TObject::kOverwrite);
    NConst2_large.at(i)->Write("",TObject::kOverwrite);
    N_small.at(i)->Write("",TObject::kOverwrite);
    N_bsmall.at(i)->Write("",TObject::kOverwrite);
    N_large.at(i)->Write("",TObject::kOverwrite);

    PTF1_large.at(i)->Write("",TObject::kOverwrite);
    PTF2_large.at(i)->Write("",TObject::kOverwrite);
    PTF1_small.at(i)->Write("",TObject::kOverwrite);
    PTF2_small.at(i)->Write("",TObject::kOverwrite);
    PTF1_bsmall.at(i)->Write("",TObject::kOverwrite);

    ECF21_large.at(i)->Write("",TObject::kOverwrite);
    ECF22_large.at(i)->Write("",TObject::kOverwrite);
    ECF21_small.at(i)->Write("",TObject::kOverwrite);
    ECF22_small.at(i)->Write("",TObject::kOverwrite);

    Lep2PT.at(i)->Write("",TObject::kOverwrite);
    LepDeltaR.at(i)->Write("",TObject::kOverwrite);

    IMsmall_DeltaRs.at(i)->Write("",TObject::kOverwrite);
    IMlarge_DeltaRl.at(i)->Write("",TObject::kOverwrite);

    PTF1s_PTF2s.at(i)->Write("",TObject::kOverwrite);
    PTF1s_PTF1l.at(i)->Write("",TObject::kOverwrite);
    PTF1bs_PTF1l.at(i)->Write("",TObject::kOverwrite);
  }

}// end BuildHistos
