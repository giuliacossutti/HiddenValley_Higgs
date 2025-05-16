/*
Author: Giulia Cossutti

Macro to fill histograms for signal and backgrounds

From inside the /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes directory run with
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250515_Task19/BuildHistos.C'("/gfsvol01/atlas/giuliac/condor/20250513_Task19/50K_signal_EJ.root", "signal_EJ")' &
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250515_Task19/BuildHistos.C'("/gfsvol01/atlas/giuliac/condor/20250513_Task19/50K_bkg_Zjets_q.root", "bkg_Zjets_q")' &
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250515_Task19/BuildHistos.C'("/gfsvol01/atlas/giuliac/condor/20250513_Task19/50K_bkg_Zjets_g.root", "bkg_Zjets_g")' &
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250515_Task19/BuildHistos.C'("/gfsvol01/atlas/giuliac/condor/20250513_Task19/50K_bkg_ZZ.root", "bkg_ZZ")' &
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250515_Task19/BuildHistos.C'("/gfsvol01/atlas/giuliac/condor/20250513_Task19/50K_bkg_HZ_SM.root", "bkg_HZ_SM")' &

Pre-selection on events:
- At least 2 small jets
- 2 charged leptons with same flavour, opposite charge, | Invmass - mZ | < tolerance = 10 GeV, ZpT > 20 GeV
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
// Function to calculate 2-point Energy Correlation Function of a jet
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
//----------- END FUNCTIONS ---------------------


void BuildHistos(const char *inFile, const char *outFile)
{

  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inFile);

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
  

  // Book histograms
  TH1F *Invmass_2small = new TH1F("Invariant_mass_of_leading_and_subleading_small_jets", "", 100, 0.0, 400.0);
  TH1F *Invmass_2large = new TH1F("Invariant_mass_of_leading_and_subleading_large_jets", "", 100, 0.0, 400.0);
  TH1F *DeltaR_small = new TH1F("DeltaR_between_leading_and_subleading_small_jets", "", 100, 0., 10.);
  TH1F *DeltaPhi_small = new TH1F("DeltaPhi_between_leading_and_subleading_small_jets", "", 110, -1., 10.);
  TH1F *JetPT1_small = new TH1F("leading_small_jet_pt", "", 75, 0.0, 150.0);
  TH1F *JetPT2_small = new TH1F("subleading_small_jet_pt", "", 75, 0.0, 150.0);
  TH1F *JetPT1_large = new TH1F("leading_large_jet_pt", "", 75, 0.0, 150.0);
  TH1F *JetPT2_large = new TH1F("subleading_large_jet_pt", "", 75, 0.0, 150.0);
  TH1F *NConst1_large = new TH1F("number_of_constituents_of_leading_large_jet", "", 30, -0.5, 29.5);
  TH1F *NConst2_large = new TH1F("number_of_constituents_of_subleading_large_jet", "", 30, -0.5, 29.5);
  TH1F *N_small = new TH1F("number_of_small_jets", "", 30, -0.5, 29.5);
  TH1F *N_bsmall = new TH1F("number_of_btagged_small_jets", "", 30, -0.5, 29.5);
  TH1F *N_large = new TH1F("number_of_large_jets", "", 30, -0.5, 29.5);
  TH1F *PTF1_large = new TH1F("Prompt_Track_Fraction_of_leading_large_jet", "", 100, 0., 1.);
  TH1F *PTF2_large = new TH1F("Prompt_Track_Fraction_of_subleading_large_jet", "", 100, 0., 1.);
  TH1F *PTF1_small = new TH1F("Prompt_Track_Fraction_of_leading_small_jet", "", 100, 0., 1.);
  TH1F *PTF2_small = new TH1F("Prompt_Track_Fraction_of_subleading_small_jet", "", 100, 0., 1.);
  TH1F *PTF1_bsmall = new TH1F("Prompt_Track_Fraction_of_leading_btagged_small_jet", "", 100, 0., 1.);
  TH1F *ECF21_large = new TH1F("2point_Energy_Correlation_div_pT_of_leading_large_jet", "", 100, 0., 200.);
  TH1F *ECF22_large = new TH1F("2point_Energy_Correlation_div_pT_of_subleading_large_jet", "", 100, 0., 200.);

  TH2F *IMsmall_IMlarge = new TH2F("IMsmall_vs_IMlarge", "", 100, 0.0, 400.0, 100, 0.0, 400.0);
  TH2F *JetPT1s_JetPT2s = new TH2F("JetPT1small_vs_JetPT2small", "", 75, 0.0, 150.0, 75, 0.0, 150.0);
  TH2F *JetPT1l_JetPT2l = new TH2F("JetPT1large_vs_JetPT2large", "", 75, 0.0, 150.0, 75, 0.0, 150.0);
  TH2F *DeltaRs_DeltaPhis = new TH2F("DeltaRsmall_vs_DeltaPhismall", "", 100, 0., 10., 260, -1., 10.);

  // Counter of used events
  Int_t NofEv = 0;

  // Tolerance on Invariant mass near Z mass
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

     // Counter of small
     Int_t n_smalljets = 0;
     
     // Loop over all small jets
     for(Int_t jetentry = 0; jetentry < branchJet->GetEntries(); ++jetentry)
     {
       Jet *jet = (Jet*) branchJet->At(jetentry);
       //if( abs(jet->Eta) >= 2.5 ) continue;
       n_smalljets += 1;
     }//end loop over all small jets

     if( (n_smalljets < 2) ) continue;



     // Preselection on charged leptons from Z decay

     // Is there acceptable Z -> l+ l- ?
     Bool_t isZ2l = kFALSE; 

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
             }
           }
         }
       }
     }

     if( isZ2l == kFALSE ) continue;

     // Increment counter of events passing the preselection
     NofEv += 1;
     
    //-------------------- End of Preselection -------------------------

     // jet counter
     Int_t index_small = 0;
     Int_t index_bsmall = 0;
     Int_t index_large = 0;

     Int_t index_small1 = 0;
     Int_t index_bsmall1 = 0;
     Int_t index_small2 = 0;

     Int_t index_large1 = 0;
     Int_t index_large2 = 0;


     // Loop over all small jets
     for(Int_t jetentry = 0; jetentry < branchJet->GetEntries(); ++jetentry)
     {
       Jet *jet = (Jet*) branchJet->At(jetentry);
       //if( abs(jet->Eta) >= 2.5 ) continue;
       index_small += 1;
       if( jet->BTag == 1 ) index_bsmall += 1;

       // Indexes of leading and subleading small jets
       if(index_small == 1){
         index_small1 = jetentry;
       }else if(index_small == 2){
         index_small2 = jetentry;
       }
       if(index_bsmall == 1){
         index_bsmall1 = jetentry;
       }
     }//end loop over all small jets

     // Loop over all large jets
     for(Int_t jetentry = 0; jetentry < branchReclusteredJet->GetEntries(); ++jetentry)
     {
       Jet *jet = (Jet*) branchReclusteredJet->At(jetentry);
       //if( abs(jet->Eta) >= 1.8 ) continue;
       index_large += 1;

       // Indexes of leading and subleading large jets
       if(index_large == 1){
         index_large1 = jetentry;
       }else if(index_large == 2){
         index_large2 = jetentry;
       }
     }//end loop over all large jets



     // Fill histograms

     N_small->Fill(index_small);
     N_bsmall->Fill(index_bsmall);
     N_large->Fill(index_large);

     if(index_bsmall >= 1){
       Jet *smallbjet1 = (Jet*) branchJet->At(index_bsmall1);
       PTF1_bsmall->Fill( PTFsmall(treeReader, entry, branchTrack, smallbjet1, zPV) );
     }

     if(index_small >= 2){
       Jet *smalljet1 = (Jet*) branchJet->At(index_small1);
       Jet *smalljet2 = (Jet*) branchJet->At(index_small2);

       TLorentzVector *smalljet1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       smalljet1P4->SetPtEtaPhiM(smalljet1->PT,smalljet1->Eta,smalljet1->Phi,smalljet1->Mass);
       TLorentzVector *smalljet2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       smalljet2P4->SetPtEtaPhiM(smalljet2->PT,smalljet2->Eta,smalljet2->Phi,smalljet2->Mass);

       //Fill 1D histograms
       Invmass_2small->Fill( (*smalljet1P4 + *smalljet2P4).M()  );

       DeltaR_small->Fill( DeltaR(smalljet1, smalljet2) );
       DeltaPhi_small->Fill( abs(DeltaPhi(smalljet1, smalljet2)) );

       JetPT1_small->Fill( smalljet1->PT );
       JetPT2_small->Fill( smalljet2->PT );

       PTF1_small->Fill( PTFsmall(treeReader, entry, branchTrack, smalljet1, zPV) );
       PTF2_small->Fill( PTFsmall(treeReader, entry, branchTrack, smalljet2, zPV) );

       //Fill 2D histograms
       JetPT1s_JetPT2s->Fill( smalljet1->PT , smalljet2->PT );
       DeltaRs_DeltaPhis->Fill( DeltaR(smalljet1, smalljet2) , abs(DeltaPhi(smalljet1, smalljet2)) );
     }

     if(index_large >= 2){
       Jet *largejet1 = (Jet*) branchReclusteredJet->At(index_large1);
       Jet *largejet2 = (Jet*) branchReclusteredJet->At(index_large2);

       TLorentzVector *largejet1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       largejet1P4->SetPtEtaPhiM(largejet1->PT,largejet1->Eta,largejet1->Phi,largejet1->Mass);
       TLorentzVector *largejet2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       largejet2P4->SetPtEtaPhiM(largejet2->PT,largejet2->Eta,largejet2->Phi,largejet2->Mass);

       //Fill 1D histograms
       Invmass_2large->Fill( (*largejet1P4 + *largejet2P4).M()  );

       JetPT1_large->Fill( largejet1->PT );
       JetPT2_large->Fill( largejet2->PT );

       NConst1_large->Fill(largejet1->Constituents.GetEntriesFast());       
       NConst2_large->Fill(largejet2->Constituents.GetEntriesFast());       

       PTF1_large->Fill( PTF(treeReader, entry, branchTrack, largejet1, zPV) );
       PTF2_large->Fill( PTF(treeReader, entry, branchTrack, largejet2, zPV) );

       ECF21_large->Fill( ECF2(treeReader, entry, branchTrack, largejet1) / largejet1->PT );
       ECF22_large->Fill( ECF2(treeReader, entry, branchTrack, largejet2) / largejet2->PT );

       //Fill 2D histograms
       JetPT1l_JetPT2l->Fill( largejet1->PT , largejet2->PT );

       if(index_small >= 2){
         Jet *smalljet1 = (Jet*) branchJet->At(index_small1);
         Jet *smalljet2 = (Jet*) branchJet->At(index_small2);

         TLorentzVector *smalljet1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
         smalljet1P4->SetPtEtaPhiM(smalljet1->PT,smalljet1->Eta,smalljet1->Phi,smalljet1->Mass);
         TLorentzVector *smalljet2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
         smalljet2P4->SetPtEtaPhiM(smalljet2->PT,smalljet2->Eta,smalljet2->Phi,smalljet2->Mass);

         IMsmall_IMlarge->Fill(  (*smalljet1P4 + *smalljet2P4).M() , (*largejet1P4 + *largejet2P4).M() );
       }

     }

  }// end loop over all events

//------------------------------------------------

  //cout << "Number of used Events: " << NofEv << endl;

//-------------------------------------------------------
  // write all plots on file
  TString outname;
  outname.Form("/gfsvol01/atlas/giuliac/plots_and_outputs/20250515_Task19/%s.root",outFile);

  TFile *f = new TFile(outname,"RECREATE");

  Invmass_2small->Write("",TObject::kOverwrite);
  Invmass_2large->Write("",TObject::kOverwrite);
  DeltaR_small->Write("",TObject::kOverwrite);
  DeltaPhi_small->Write("",TObject::kOverwrite);
  JetPT1_small->Write("",TObject::kOverwrite);
  JetPT2_small->Write("",TObject::kOverwrite);
  JetPT1_large->Write("",TObject::kOverwrite);
  JetPT2_large->Write("",TObject::kOverwrite);
  NConst1_large->Write("",TObject::kOverwrite);
  NConst2_large->Write("",TObject::kOverwrite);
  N_small->Write("",TObject::kOverwrite);
  N_bsmall->Write("",TObject::kOverwrite);
  N_large->Write("",TObject::kOverwrite);
  PTF1_large->Write("",TObject::kOverwrite);
  PTF2_large->Write("",TObject::kOverwrite);
  PTF1_small->Write("",TObject::kOverwrite);
  PTF2_small->Write("",TObject::kOverwrite);
  PTF1_bsmall->Write("",TObject::kOverwrite);
  ECF21_large->Write("",TObject::kOverwrite);
  ECF22_large->Write("",TObject::kOverwrite);

  DeltaRs_DeltaPhis->Write("",TObject::kOverwrite);
  JetPT1s_JetPT2s->Write("",TObject::kOverwrite);
  JetPT1l_JetPT2l->Write("",TObject::kOverwrite);
  IMsmall_IMlarge->Write("",TObject::kOverwrite);

}// end BuildHistos
