/*
Author: Giulia Cossutti

Macro to plot number of jets with R=0.4 and reclustered jets with R=1.0
Macro to plot Invariant Mass of jets with R=0.4 and reclustered jets with R=1.0
Macro to plot DeltaR, DeltaEta, DeltaPhi between jets with R=0.4 and reclustered jets with R=1.0
Macro to plot Jet PT (various combinations)

From inside the /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes directory run with
root -l /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250520_Task20/Reclustering.C'("/gfsvol01/atlas/giuliac/condor/20250520_Task20/50K_signal_EJ.root","signal_EJ",kTRUE)'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "/gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/classes/DelphesClasses.h"
#include "/gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/external/ExRootAnalysis/ExRootTreeReader.h"
#include "/gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/external/ExRootAnalysis/ExRootResult.h"
#endif

//------------------------------------------------------------------------------

//---------- FUNCTIONS -------------------
//===================================
// Function to customize histograms:
// - Scale by number of entries and bin width
// - Customize markers and lines
void CustomizeHist(TH1F* hist, Style_t style, Color_t color, Bool_t scale = kFALSE){
  hist->Sumw2();
  if(scale == kTRUE) hist->Scale(1.0/(hist->GetEntries() * hist->GetXaxis()->GetBinWidth(2)));
  hist->SetMarkerStyle(style);
  hist-> SetMarkerColor(color);
  hist-> SetLineColor(color);
  hist->SetMarkerSize(0.8);
  hist->GetXaxis()->SetTitleSize(.045);
  hist->GetYaxis()->SetTitleSize(.045);
}
//=================================
// Function to customize legends:
// - Customize text size of header and body
// - No border
void CustomizeLeg(TLegend* legend){
  TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
  header->SetTextSize(0.03);
  gStyle->SetLegendTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetEntrySeparation(0.01);
}
//===============================
// Function to compute DeltaR between 2 jets
Double_t DeltaR(Jet* jet1, Jet* jet2){
 return sqrt( pow(jet1->Eta - jet2->Eta,2) + pow(jet1->Phi - jet2->Phi,2) );
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
// Function to loop over constituents of a jet
void LoopOverConstituents(Jet* jet, TString string){
         for(Int_t j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
          {
            // Define objects for loop on constituents
            TObject *object;
            GenParticle *particle;
            //Electron *electron;
            //Photon *photon;
            Muon *muon;
            Track *track;
            Tower *tower;
            Jet *jet_in;

            object = jet->Constituents.At(j);

            // Check if the constituent is accessible
            if(object == 0) continue;
  
            cout << string << " JET  " <<jet->Constituents.GetEntriesFast() << " constituents. Constituent is a " << object->IsA() << endl;

            if(object->IsA() == GenParticle::Class())
            {
              particle = (GenParticle*) object;
               cout << "    GenPart pt: " << particle->PT << ", eta: " << particle->Eta << ", phi: " << particle->Phi << endl;
            }
            else if(object->IsA() == Track::Class())
            {
              track = (Track*) object;
               cout << "    Track pt: " << track->PT << ", eta: " << track->Eta << ", phi: " << track->Phi << endl;
            }
            else if(object->IsA() == Tower::Class())
            {
              tower = (Tower*) object;
               cout << "    Tower pt: " << tower->ET << ", eta: " << tower->Eta << ", phi: " << tower->Phi << endl;
            }
            else if(object->IsA() == Jet::Class())
            {
              jet_in = (Jet*) object;
               cout << "    Jet pt: " << jet_in->PT << ", eta: " << jet_in->Eta << ", phi: " << jet_in->Phi << endl;
            }
            else if(object->IsA() == Muon::Class())
            {
              muon = (Muon*) object;
               cout << "    Muon pt: " << muon->PT  << endl;
            }
          }// end loop over constituents of large jets
}// end function to loop over constituents of a jet
//===============================
//---------- END FUNCTIONS -------------------


void Reclustering(const char *inputFile, const char * outFile,Bool_t Norm=kFALSE)
{
  // No Event selection or Small and Large Jet definition: 

  // Output name:
  TString outfolder;
  outfolder.Form("/gfsvol01/atlas/giuliac/plots_and_outputs/20250520_Task20/%s_",outFile);

  // Normalisation of histograms
  Bool_t NjetsNorm = Norm;
  Bool_t InvMassNorm = Norm;
  Bool_t DeltaRNorm = Norm;
  Bool_t DeltaEtaNorm = Norm;
  Bool_t DeltaPhiNorm = Norm;
  Bool_t JetPTNorm = Norm;

  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

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
  TH1F *Njets_small = new TH1F("number_of_jets_with_R0.4", "", 30, -0.5, 29.5);
  TH1F *Nbjets_small = new TH1F("number_of_btagged_jets_with_R0.4", "", 30, -0.5, 29.5);
  TH1F *Njets_large = new TH1F("number_of_reclustered_jets_with_R1.0", "", 30, -0.5, 29.5);
  TH1F *Njets_in_large = new TH1F("number_of_constituents_of_leading_and_subleading_jets_with_R1.0", "", 30, -0.5, 29.5);
  TH1F *Njets_in_small = new TH1F("number_of_constituents_of_leading_and_subleading_jets_with_R0.4", "", 30, -0.5, 29.5);

  TH1F *Invmass_2small = new TH1F("Invariant_mass_of_Leading_and_Subleading_small_jets", "", 100, 0.0, 1500.0);
  TH1F *Invmass_4small = new TH1F("Invariant_mass_of_first four_small_jets", "", 100, 0.0, 1500.0);
  TH1F *Invmass_2large = new TH1F("Invariant_mass_of_Leading_and_Subleading_large_jets", "", 100, 0.0, 1500.0);

  TH1F *DeltaR_large = new TH1F("DeltaR_between_leading_and_subleading_jets_with_R1.0", "", 100, 0., 10.);
  TH1F *DeltaR_large2 = new TH1F("DeltaR_between_leading_and_subleading_jets_with_R1.0_if_exactly_2_jets", "", 100, 0., 10.);
  TH1F *DeltaR_small = new TH1F("DeltaR_between_leading_and_subleading_jets_with_R0.4", "", 100, 0., 10.);

  TH1F *DeltaEta_large = new TH1F("DeltaEta_between_leading_and_subleading_jets_with_R1.0", "", 220, -8., 14.);
  TH1F *DeltaEta_large2 = new TH1F("DeltaEta_between_leading_and_subleading_jets_with_R1.0_if_exactly_2_jets", "", 220, -8., 14.);
  TH1F *DeltaEta_small = new TH1F("DeltaEta_between_leading_and_subleading_jets_with_R0.4", "", 220, -8., 14.);

  TH1F *DeltaPhi_large = new TH1F("DeltaPhi_between_leading_and_subleading_jets_with_R1.0", "", 240, -8., 16.);
  TH1F *DeltaPhi_large2 = new TH1F("DeltaPhi_between_leading_and_subleading_jets_with_R1.0_if_exactly_2_jets", "", 240, -8., 16.);
  TH1F *DeltaPhi_small = new TH1F("DeltaPhi_between_leading_and_subleading_jets_with_R0.4", "", 240, -8., 16.);

  Double_t PTend = 800.;
  TH1F *JetPT_small = new TH1F("all_small_jet_pt", "", PTend/10., 0.0, PTend);
  TH1F *JetPT1_small = new TH1F("leading_small_jet_pt", "", PTend/10., 0.0, PTend);
  TH1F *JetPT2_small = new TH1F("subleading_small_jet_pt", "", PTend/10., 0.0, PTend);
  TH1F *JetPT12_small = new TH1F("pt_of_first_two_small_jets", "", PTend/10., 0.0, PTend);
  TH1F *JetPT4_small = new TH1F("pt_of_first_four_small_jets", "", PTend/10., 0.0, PTend);

  TH1F *JetPT_large = new TH1F("all_large_jet_pt", "", PTend/10., 0.0, PTend);
  TH1F *JetPT1_large = new TH1F("leading_large_jet_pt", "", PTend/10., 0.0, PTend);
  TH1F *JetPT2_large = new TH1F("subleading_large_jet_pt", "", PTend/10., 0.0, PTend);
  TH1F *JetPT12_large = new TH1F("pt_of_first_two_large_jets", "", PTend/10., 0.0, PTend);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

     // jet counter
     Int_t index_small = 0;
     Int_t index_large = 0;

     // jetentry of leading and subleading jet
     Int_t index_small1 = 0;
     Int_t index_small2 = 0;
     Int_t index_small3 = 0;
     Int_t index_small4 = 0;

     Int_t index_large1 = 0;
     Int_t index_large2 = 0;

     Int_t n_smallbjets = 0;

     // Loop over all small jets
     for(Int_t jetentry = 0; jetentry < branchJet->GetEntries(); ++jetentry)
     {
       Jet *jet = (Jet*) branchJet->At(jetentry);
       //if( abs(jet->Eta) >= 2.5 ) continue;
       index_small += 1;
       if( jet->BTag == 1 ) n_smallbjets += 1;

       // Fill all small jet pT histogram
       JetPT_small->Fill( jet->PT );

       // Pt of leading or subleading small jet
       if(index_small == 1) JetPT1_small->Fill( jet->PT );
       if(index_small == 2) JetPT2_small->Fill( jet->PT );

       // Indexes of leading and subleading jets
       if(index_small == 1) index_small1 = jetentry;
       if(index_small == 2) index_small2 = jetentry;
       if(index_small == 3) index_small3 = jetentry;
       if(index_small == 4) index_small4 = jetentry;

       // Loop over leading and subleading small jet constituents
       if( index_small == 1 || index_small == 2 ){
         //if(entry <= 5) LoopOverConstituents(jet, "SMALL");   
         Njets_in_small->Fill(jet->Constituents.GetEntriesFast());     
       }//end if

     }//end loop over all small jets

     // Loop over all large jets
     for(Int_t jetentry = 0; jetentry < branchReclusteredJet->GetEntries(); ++jetentry)
     {
       Jet *jet = (Jet*) branchReclusteredJet->At(jetentry);
       //if( abs(jet->Eta) >= 1.8 ) continue;
       index_large += 1;

       // Fill all large jet pT histogram
       JetPT_large->Fill( jet->PT );

       // Pt of leading or subleading large jet
       if(index_large == 1) JetPT1_large->Fill( jet->PT );
       if(index_large == 2) JetPT2_large->Fill( jet->PT );

       // Indexes of leading and subleading jets
       if(index_large == 1) index_large1 = jetentry;
       if(index_large == 2) index_large2 = jetentry;
       
       // Loop over leading and subleading large jet constituents
       if (index_large == 1 || index_large == 2){   
         //if(entry <= 5) LoopOverConstituents(jet, "LARGE");   
         Njets_in_large->Fill(jet->Constituents.GetEntriesFast());
       }//end if

     }//end loop over all large jets

     // Fill number of jets histograms
     Njets_small->Fill(index_small);
     Nbjets_small->Fill(n_smallbjets);
     Njets_large->Fill(index_large);

     // Invariant masses, DeltaR, DeltaPhi, DeltaEta
     if(index_small >= 2){
       Jet *smalljet1 = (Jet*) branchJet->At(index_small1);
       Jet *smalljet2 = (Jet*) branchJet->At(index_small2);

       TLorentzVector *smalljet1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       smalljet1P4->SetPtEtaPhiM(smalljet1->PT,smalljet1->Eta,smalljet1->Phi,smalljet1->Mass);
       TLorentzVector *smalljet2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       smalljet2P4->SetPtEtaPhiM(smalljet2->PT,smalljet2->Eta,smalljet2->Phi,smalljet2->Mass);

       Invmass_2small->Fill( (*smalljet1P4 + *smalljet2P4).M()  );
       DeltaR_small->Fill( DeltaR(smalljet1, smalljet2) );
       DeltaEta_small->Fill( DeltaEta(smalljet1, smalljet2) );
       DeltaPhi_small->Fill( DeltaPhi(smalljet1, smalljet2) );

       // Pt of sum of first two small jets
       JetPT12_small->Fill( (*smalljet1P4 + *smalljet2P4).Pt()  );

       if(index_small >= 4){
         Jet *smalljet3 = (Jet*) branchJet->At(index_small3);
         Jet *smalljet4 = (Jet*) branchJet->At(index_small4);

         TLorentzVector *smalljet3P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
         smalljet3P4->SetPtEtaPhiM(smalljet3->PT,smalljet3->Eta,smalljet3->Phi,smalljet3->Mass);
         TLorentzVector *smalljet4P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
         smalljet4P4->SetPtEtaPhiM(smalljet4->PT,smalljet4->Eta,smalljet4->Phi,smalljet4->Mass);

         Invmass_4small->Fill( (*smalljet1P4 + *smalljet2P4 + *smalljet3P4 + *smalljet4P4).M()  )
;
         // Pt of sum of first two small jets
         JetPT4_small->Fill( (*smalljet1P4 + *smalljet2P4 + *smalljet3P4 + *smalljet4P4).Pt()  );
       }
     }

     if(index_large >= 2){
       Jet *largejet1 = (Jet*) branchReclusteredJet->At(index_large1);
       Jet *largejet2 = (Jet*) branchReclusteredJet->At(index_large2);

       TLorentzVector *largejet1P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       largejet1P4->SetPtEtaPhiM(largejet1->PT,largejet1->Eta,largejet1->Phi,largejet1->Mass);
       TLorentzVector *largejet2P4 = new TLorentzVector(1.0,1.0,1.0,1.0);
       largejet2P4->SetPtEtaPhiM(largejet2->PT,largejet2->Eta,largejet2->Phi,largejet2->Mass);

       Invmass_2large->Fill( (*largejet1P4 + *largejet2P4).M()  );
       DeltaR_large->Fill( DeltaR(largejet1, largejet2) );
       DeltaEta_large->Fill( DeltaEta(largejet1, largejet2) );
       DeltaPhi_large->Fill( DeltaPhi(largejet1, largejet2) );
       if(index_large == 2) DeltaR_large2->Fill( DeltaR(largejet1, largejet2) );
       if(index_large == 2) DeltaEta_large2->Fill( DeltaEta(largejet1, largejet2) );
       if(index_large == 2) DeltaPhi_large2->Fill( DeltaPhi(largejet1, largejet2) );

       // Pt of sum of first two large jets
       JetPT12_large->Fill( (*largejet1P4 + *largejet2P4).Pt()  );
     }


  }//end loop over all events

//-----------------------------------------------------------------

  // Canva to draw Number of jets histograms
  TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);

  //Bool_t NjetsNorm = kTRUE;
  //Bool_t NjetsNorm = kFALSE;
  CustomizeHist(Njets_small, kFullSquare, kBlue, NjetsNorm);
  CustomizeHist(Nbjets_small, kFullSquare, kViolet, NjetsNorm);
  CustomizeHist(Njets_large, kFullCircle, kRed, NjetsNorm);
  CustomizeHist(Njets_in_large, kFullCircle, kGreen, NjetsNorm);
  CustomizeHist(Njets_in_small, kFullSquare, kBlack, NjetsNorm);

  Njets_large->GetXaxis()->SetTitle("Number of jets or constituents");
  if(NjetsNorm == kTRUE){Njets_large->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{dN_{jets or constituents}}");}
  else{Njets_large->GetYaxis()->SetTitle("a.u.");}

  Njets_in_large->GetXaxis()->SetTitle("Number of jets or constituents");
  if(NjetsNorm == kTRUE){Njets_in_large->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{dN_{jets or constituents}}");}
  else{Njets_in_large->GetYaxis()->SetTitle("a.u.");}

  Nbjets_small->GetXaxis()->SetTitle("Number of jets or constituents");
  if(NjetsNorm == kTRUE){Nbjets_small->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{dN_{jets or constituents}}");}
  else{Nbjets_small->GetYaxis()->SetTitle("a.u.");}

  // Show resulting histograms
  
  if( (strcmp(outFile,"ReclpTmin40")==0) || (strcmp(outFile,"bbbar")==0) || (strcmp(outFile,"signal_EJ")==0) ){
    Njets_in_large->Draw("E,HIST");
    Njets_large->Draw("E,HIST same");
    Njets_small->Draw("E,HIST same");
    Njets_in_small->Draw("E,HIST same");
    Nbjets_small->Draw("E,HIST same");
  }else{
    Njets_large->Draw("E,HIST");
    Njets_small->Draw("E,HIST same");
    Njets_in_large->Draw("E,HIST same");
    Njets_in_small->Draw("E,HIST same");
    Nbjets_small->Draw("E,HIST same");
  }

  // Legend
  TLegend* legend = new TLegend(0.4599,0.385417,0.807018,0.864583);

  legend->SetHeader(TString::Format("#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, %s}{#splitline{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV}{Number of Jets or Constituents with Reclustering}}",outFile),"C");

  legend->AddEntry(Njets_small,"Small (R=0.4)");
  legend->AddEntry(Nbjets_small,"Small (R=0.4) b-tagged");
  legend->AddEntry(Njets_large,"Large (R=1.0)");
  legend->AddEntry(Njets_in_small,"Constituents of L or Subl Small");
  legend->AddEntry(Njets_in_large,"Constituents of L or Subl Large");
  CustomizeLeg(legend);
  legend->Draw();

  if(NjetsNorm==kFALSE){
    c1->SaveAs(outfolder + "Njets.pdf"); 
  }else if(NjetsNorm==kTRUE){
    c1->SaveAs(outfolder + "Njets_Norm.pdf"); 
  }
  
  //c1->Close();

//-------------------------------------------------------

  // Canva to draw InvMass histograms
  TCanvas *c2 = new TCanvas("c2", "c2", 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();
  //gPad->SetLogx();

  //Bool_t InvMassNorm = kTRUE;
  //Bool_t InvMassNorm = kFALSE;
  CustomizeHist(Invmass_2small, kFullSquare, kBlue, InvMassNorm);
  CustomizeHist(Invmass_4small, kFullSquare, kBlack, InvMassNorm);
  CustomizeHist(Invmass_2large, kFullCircle, kRed, InvMassNorm);

  Invmass_2small->GetXaxis()->SetTitle("Invariant Mass [GeV/c^2]");
  if(InvMassNorm == kTRUE){Invmass_2small->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{dM_{inv}} [c^2/GeV]");}
  else{Invmass_2small->GetYaxis()->SetTitle("a.u.");}

  // Show resulting histograms
  Invmass_2small->Draw("E,HIST");
  Invmass_4small->Draw("E,HIST same");
  Invmass_2large->Draw("E,HIST same");

  // Legend
  TLegend* legendIM = new TLegend(0.547619,0.5625,0.860902,0.875);

  legendIM->SetHeader(TString::Format("#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, %s}{#splitline{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV}{Invariant Mass of Jets}}",outFile),"C");

  legendIM->AddEntry(Invmass_2small,"First 2 Small (R=0.4)");
  legendIM->AddEntry(Invmass_4small,"First 4 Small (R=0.4)");
  legendIM->AddEntry(Invmass_2large,"First 2 Large (R=1.0)");
  CustomizeLeg(legendIM);
  legendIM->Draw();

  if(InvMassNorm==kFALSE){
    c2->SaveAs(outfolder + "InvMass.pdf"); 
  }else if(InvMassNorm==kTRUE){
    c2->SaveAs(outfolder + "InvMass_Norm.pdf"); 
  }
  
  //c2->Close();

//----------------------------------------------

  // Canva to draw DeltaR between jets histograms
  TCanvas *c3 = new TCanvas("c3", "c3", 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();

  //Bool_t DeltaRNorm = kTRUE;
  //Bool_t DeltaRNorm = kFALSE;
  CustomizeHist(DeltaR_large, kFullSquare, kBlue, DeltaRNorm);
  CustomizeHist(DeltaR_large2, kFullSquare, kBlack, DeltaRNorm);
  CustomizeHist(DeltaR_small, kFullCircle, kRed, DeltaRNorm);

  DeltaR_large->GetXaxis()->SetTitle("#Delta R");
  if(DeltaRNorm == kTRUE){DeltaR_large->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d #Delta R}");}
  else{DeltaR_large->GetYaxis()->SetTitle("a.u.");}

  DeltaR_small->GetXaxis()->SetTitle("#Delta R");
  if(DeltaRNorm == kTRUE){DeltaR_small->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d #Delta R}");}
  else{DeltaR_small->GetYaxis()->SetTitle("a.u.");}

  // Show resulting histograms
  if( Norm==kTRUE ){
    DeltaR_large->Draw("E,HIST");
    DeltaR_small->Draw("E,HIST same");
    DeltaR_large2->Draw("E,HIST same");
  }else{
    DeltaR_small->Draw("E,HIST");
    DeltaR_large->Draw("E,HIST same");
    DeltaR_large2->Draw("E,HIST same");
  }

  // Legend
  TLegend* legendDR = new TLegend(0.588972,0.567708,0.874687,0.876736);

  legendDR->SetHeader(TString::Format("#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, %s}{#splitline{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV}{#Delta R between Jets}}",outFile),"C");

  legendDR->AddEntry(DeltaR_small,"First 2 Small (R=0.4)");
  legendDR->AddEntry(DeltaR_large,"First 2 Large (R=1.0)");
  legendDR->AddEntry(DeltaR_large2,"Only 2 Large (R=1.0)");
  CustomizeLeg(legendDR);
  legendDR->Draw();

  if(DeltaRNorm==kFALSE){
    c3->SaveAs(outfolder + "DeltaR.pdf"); 
  }else if(DeltaRNorm==kTRUE){
    c3->SaveAs(outfolder + "DeltaR_Norm.pdf"); 
  }

  //c3->Close();

//----------------------------------------------

  // Canva to draw DeltaEta between jets histograms
  TCanvas *c4 = new TCanvas("c4", "c4", 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();

  //Bool_t DeltaEtaNorm = kTRUE;
  //Bool_t DeltaEtaNorm = kFALSE;
  CustomizeHist(DeltaEta_large, kFullSquare, kBlue, DeltaEtaNorm);
  CustomizeHist(DeltaEta_large2, kFullSquare, kBlack, DeltaEtaNorm);
  CustomizeHist(DeltaEta_small, kFullCircle, kRed, DeltaEtaNorm);

  DeltaEta_small->GetXaxis()->SetTitle("#Delta #eta");
  if(DeltaEtaNorm == kTRUE){DeltaEta_small->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d #Delta #eta}");}
  else{DeltaEta_small->GetYaxis()->SetTitle("a.u.");}

  // Show resulting histograms
  DeltaEta_small->Draw("E,HIST");
  DeltaEta_large2->Draw("E,HIST same");
  DeltaEta_large->Draw("E,HIST same");

  // Legend
  TLegend* legendDEta = new TLegend(0.56391,0.560764,0.877193,0.875);

  legendDEta->SetHeader(TString::Format("#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, %s}{#splitline{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV}{#Delta #eta between Jets}}",outFile),"C");

  legendDEta->AddEntry(DeltaEta_small,"First 2 Small (R=0.4)");
  legendDEta->AddEntry(DeltaEta_large,"First 2 Large (R=1.0)");
  legendDEta->AddEntry(DeltaEta_large2,"Only 2 Large (R=1.0)");
  CustomizeLeg(legendDEta);
  legendDEta->Draw();

  if(DeltaEtaNorm==kFALSE){
    c4->SaveAs(outfolder + "DeltaEta.pdf"); 
  }else if(DeltaEtaNorm==kTRUE){
    c4->SaveAs(outfolder + "DeltaEta_Norm.pdf"); 
  }
  
  //c4->Close();
//----------------------------------------------

  // Canva to draw DeltaPhi between jets histograms
  TCanvas *c5 = new TCanvas("c5", "c5", 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();

  //Bool_t DeltaPhiNorm = kTRUE;
  //Bool_t DeltaPhiNorm = kFALSE;
  CustomizeHist(DeltaPhi_large, kFullSquare, kBlue, DeltaPhiNorm);
  CustomizeHist(DeltaPhi_large2, kFullSquare, kBlack, DeltaPhiNorm);
  CustomizeHist(DeltaPhi_small, kFullCircle, kRed, DeltaPhiNorm);

  DeltaPhi_large->GetXaxis()->SetTitle("#Delta #varphi");
  if(DeltaPhiNorm == kTRUE){DeltaPhi_large->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d #Delta #varphi}");}
  else{DeltaPhi_large->GetYaxis()->SetTitle("a.u.");}

  DeltaPhi_small->GetXaxis()->SetTitle("#Delta #varphi");
  if(DeltaPhiNorm == kTRUE){DeltaPhi_small->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d #Delta #varphi}");}
  else{DeltaPhi_small->GetYaxis()->SetTitle("a.u.");}

  // Show resulting histograms
  if( Norm==kTRUE ){
    DeltaPhi_large->Draw("E,HIST");
    DeltaPhi_small->Draw("E,HIST same");
    DeltaPhi_large2->Draw("E,HIST same");
  }else{
    DeltaPhi_small->Draw("E,HIST");
    DeltaPhi_large->Draw("E,HIST same");
    DeltaPhi_large2->Draw("E,HIST same");
  }

  // Legend
  TLegend* legendDPhi = new TLegend(0.56391,0.560764,0.877193,0.875);

  legendDPhi->SetHeader(TString::Format("#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, %s}{#splitline{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV}{#Delta #varphi between Jets}}",outFile),"C");

  legendDPhi->AddEntry(DeltaPhi_small,"First 2 Small (R=0.4)");
  legendDPhi->AddEntry(DeltaPhi_large,"First 2 Large (R=1.0)");
  legendDPhi->AddEntry(DeltaPhi_large2,"Only 2 Large (R=1.0)");
  CustomizeLeg(legendDPhi);
  legendDPhi->Draw();

  if(DeltaPhiNorm==kFALSE){
    c5->SaveAs(outfolder + "DeltaPhi.pdf"); 
  }else if(DeltaPhiNorm==kTRUE){
    c5->SaveAs(outfolder + "DeltaPhi_Norm.pdf"); 
  }

  //c5->Close();

//-------------------------------------------------------

  // Canva to draw jet PT histograms
  TCanvas *c6 = new TCanvas("c6", "c6", 0, 0, 1000, 600);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();

  //Bool_t JetPTNorm = kTRUE;
  //Bool_t JetPTNorm = kFALSE;
  CustomizeHist(JetPT_large, kFullCircle, kPink+1, JetPTNorm);
  CustomizeHist(JetPT1_large, kFullCircle, kRed, JetPTNorm);
  CustomizeHist(JetPT2_large, kFullCircle, kOrange, JetPTNorm);
  CustomizeHist(JetPT12_large, kFullCircle, kYellow-3, JetPTNorm);
  CustomizeHist(JetPT_small, kFullSquare, kGreen, JetPTNorm);
  CustomizeHist(JetPT1_small, kFullSquare, kCyan, JetPTNorm);
  CustomizeHist(JetPT2_small, kFullSquare, kBlue, JetPTNorm);
  CustomizeHist(JetPT12_small, kFullSquare, kViolet, JetPTNorm);
  CustomizeHist(JetPT4_small, kFullSquare, kBlack, JetPTNorm);

  Float_t Mrksize = 0.5;
  JetPT_small->SetMarkerSize(Mrksize);
  JetPT1_small->SetMarkerSize(Mrksize);
  JetPT2_small->SetMarkerSize(Mrksize);
  JetPT12_small->SetMarkerSize(Mrksize);
  JetPT4_small->SetMarkerSize(Mrksize);
  JetPT_large->SetMarkerSize(Mrksize);
  JetPT1_large->SetMarkerSize(Mrksize);
  JetPT2_large->SetMarkerSize(Mrksize);
  JetPT12_large->SetMarkerSize(Mrksize);

  JetPT_small->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  if(JetPTNorm == kTRUE){JetPT_small->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{dp_{T}} [c/GeV]");}
  else{JetPT_small->GetYaxis()->SetTitle("a.u.");}

  // Show resulting histograms
  JetPT_small->Draw("E,HIST");
  JetPT1_small->Draw("E,HIST same");
  JetPT2_small->Draw("E,HIST same");
  JetPT_large->Draw("E,HIST same");
  JetPT1_large->Draw("E,HIST same");
  JetPT2_large->Draw("E,HIST same");
  JetPT12_small->Draw("E,HIST same");
  JetPT4_small->Draw("E,HIST same");
  JetPT12_large->Draw("E,HIST same");

  // Legend
  TLegend* legendPT = new TLegend(0.583166,0.315972,0.896794,0.824653);

  legendPT->SetHeader(TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, %s}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV}}{pT of Jets}",outFile),"C");

  legendPT->AddEntry(JetPT_large,"All Large (R=1.0)");
  legendPT->AddEntry(JetPT1_large,"Leading Large (R=1.0)");
  legendPT->AddEntry(JetPT2_large,"Subleading Large (R=1.0)");
  legendPT->AddEntry(JetPT12_large,"L+Subl Large (R=1.0)");
  legendPT->AddEntry(JetPT_small,"All Small (R=0.4)");
  legendPT->AddEntry(JetPT1_small,"Leading Small (R=0.4)");
  legendPT->AddEntry(JetPT2_small,"Subleading Small (R=0.4)");
  legendPT->AddEntry(JetPT12_small,"L+Subl Small (R=0.4)");
  legendPT->AddEntry(JetPT4_small,"Sum of first 4 Small (R=0.4)");
  CustomizeLeg(legendPT);
  legendPT->Draw();

  if(JetPTNorm==kFALSE){
    c6->SaveAs(outfolder + "JetPT.pdf"); 
  }else if(JetPTNorm==kTRUE){
    c6->SaveAs(outfolder + "JetPT_Norm.pdf"); 
  }

  //c6->Close();

//-------------------------------------------------------
  // write all plots on file
  if(Norm==kFALSE){
    TFile *f = new TFile(outfolder + "plots.root","RECREATE");
  }else if(Norm==kTRUE){
    TFile *f = new TFile(outfolder + "plots_Norm.root","RECREATE");
  }

  Njets_small->Write("",TObject::kOverwrite);
  Njets_large->Write("",TObject::kOverwrite);
  Njets_in_large->Write("",TObject::kOverwrite);

  Invmass_2small->Write("",TObject::kOverwrite);
  Invmass_4small->Write("",TObject::kOverwrite);
  Invmass_2large->Write("",TObject::kOverwrite);

  DeltaR_large->Write("",TObject::kOverwrite);
  DeltaR_large2->Write("",TObject::kOverwrite);
  DeltaR_small->Write("",TObject::kOverwrite);

  DeltaEta_large->Write("",TObject::kOverwrite);
  DeltaEta_large2->Write("",TObject::kOverwrite);
  DeltaEta_small->Write("",TObject::kOverwrite);

  DeltaPhi_large->Write("",TObject::kOverwrite);
  DeltaPhi_large2->Write("",TObject::kOverwrite);
  DeltaPhi_small->Write("",TObject::kOverwrite);

  JetPT_small->Write("",TObject::kOverwrite);
  JetPT1_small->Write("",TObject::kOverwrite);
  JetPT2_small->Write("",TObject::kOverwrite);
  JetPT12_small->Write("",TObject::kOverwrite);
  JetPT4_small->Write("",TObject::kOverwrite);

  JetPT_large->Write("",TObject::kOverwrite);
  JetPT1_large->Write("",TObject::kOverwrite);
  JetPT2_large->Write("",TObject::kOverwrite);
  JetPT12_large->Write("",TObject::kOverwrite);

}

