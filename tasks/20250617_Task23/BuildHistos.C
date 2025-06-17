/*
Author: Giulia Cossutti

Macro to fill histograms for signal and backgrounds
2 cathegories of events:
- A: Events with 1 Large Jet
- B: Events with >=2 Large Jets

From inside the /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes directory run with
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250617_Task23/BuildHistos.C'("signal_EJ")' &
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250617_Task23/BuildHistos.C'("bkg_Zjets_q")' &
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250617_Task23/BuildHistos.C'("bkg_Zjets_g")' &
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250617_Task23/BuildHistos.C'("bkg_ZZ")' &
root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250617_Task23/BuildHistos.C'("bkg_HZ_SM")' &

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
  

  // Book histograms
  // 2 cathegories of events:
  // - A: Events with 1 Large Jet
  // - B: Events with >=2 Large Jets 

  // 1D
  TH1F *A_N_bsmall = new TH1F("A_number_of_btagged_small_jets", "", 12, -1.5, 10.5);
  TH1F *A_DeltaR_small = new TH1F("A_DeltaR_between_leading_and_subleading_small_jets", "", 100, 0., 10.);
  TH1F *A_Lep2PT = new TH1F("A_pt_of_first_two_leptons", "", 150, 0.0, 600.);

  TH1F *B_N_bsmall = (TH1F*)A_N_bsmall->Clone("B_number_of_btagged_small_jets");
  TH1F *B_DeltaR_small = (TH1F*)A_DeltaR_small->Clone("B_DeltaR_between_leading_and_subleading_small_jets");
  TH1F *B_Lep2PT = (TH1F*)A_Lep2PT->Clone("B_pt_of_first_two_leptons");

  // 2D
  TH2F *A_Nsmall_Nbsmall = new TH2F("A_Nsmall_vs_Nbsmall", "", 13, -0.5, 12.5, 12, -1.5, 10.5);

  TH1F *B_Nsmall_Nbsmall = (TH1F*)A_Nsmall_Nbsmall->Clone("B_Nsmall_vs_Nbsmall");

  // Counter of used events
  Int_t NofEv = 0;

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



     // Count number of large jets to determine cathegory of the event

     // Counter of large
     Int_t n_largejets = 0;
     
     // Loop over all large jets
     for(Int_t jetentry = 0; jetentry < branchReclusteredJet->GetEntries(); ++jetentry)
     {
       Jet *jet = (Jet*) branchReclusteredJet->At(jetentry);
       //if( abs(jet->Eta) >= 1.8 ) continue;
       n_largejets += 1;
     }//end loop over all large jets



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

             // Fill histograms about first 2 leptons when Z is tagged
             if(n_largejets == 1){
               A_Lep2PT->Fill( (*el1P4 + *el2P4).Pt() );
             }else if(n_largejets >= 2){
               B_Lep2PT->Fill( (*el1P4 + *el2P4).Pt() );
             }
 
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

               // Fill histograms about first 2 leptons when Z is tagged
               if(n_largejets == 1){
                 A_Lep2PT->Fill( (*mu1P4 + *mu2P4).Pt() );
               }else if(n_largejets >= 2){
                 B_Lep2PT->Fill( (*mu1P4 + *mu2P4).Pt() );
               }

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

     Int_t index_small1 = 0;
     Int_t index_small2 = 0;


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
     }//end loop over all small jets



     // Fill histograms

     if(n_largejets == 1){
       A_N_bsmall->Fill(index_bsmall);
       A_Nsmall_Nbsmall->Fill(index_small,index_bsmall);
     }else if(n_largejets >= 2){
       B_N_bsmall->Fill(index_bsmall);
       B_Nsmall_Nbsmall->Fill(index_small,index_bsmall);
     }

     if(index_small >= 2){
       Jet *smalljet1 = (Jet*) branchJet->At(index_small1);
       Jet *smalljet2 = (Jet*) branchJet->At(index_small2);

       //Fill 1D histograms

       if(n_largejets == 1){
         A_DeltaR_small->Fill( DeltaR(smalljet1, smalljet2) );
       }else if(n_largejets >= 2){
         B_DeltaR_small->Fill( DeltaR(smalljet1, smalljet2) );
       }
     }

  }// end loop over all events

//------------------------------------------------

  //cout << "Number of used Events for " << inFile << " : " << NofEv << endl;

//-------------------------------------------------------
  // write all plots on file
  TString outname;
  outname.Form("/gfsvol01/atlas/giuliac/plots_and_outputs/20250617_Task23/%s.root",inFile);

  TFile *f = new TFile(outname,"RECREATE");

  // 1D
  A_DeltaR_small->Write("",TObject::kOverwrite);
  B_DeltaR_small->Write("",TObject::kOverwrite);

  A_N_bsmall->Write("",TObject::kOverwrite);
  B_N_bsmall->Write("",TObject::kOverwrite);

  A_Lep2PT->Write("",TObject::kOverwrite);
  B_Lep2PT->Write("",TObject::kOverwrite);

  // 2D
  A_Nsmall_Nbsmall->Write("",TObject::kOverwrite);
  B_Nsmall_Nbsmall->Write("",TObject::kOverwrite);

}// end BuildHistos
