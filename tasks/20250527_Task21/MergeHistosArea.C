/*
Author: Giulia Cossutti

Macro to merge histograms for signal and backgrounds.
Bkgs are scaled with pythia's calculated cross-sections or measured cross-sections, 
Signal is scaled to have same area as sum of bkgs.
Histograms are saved in files.

From inside the /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes directory run with
root -l /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250527_Task21/MergeHistosArea.C 

*/

//------------------------------------------------------------------------------

//------------ FUNCTIONS ---------------------
//===================================
// Function to customize signal histogram:
void CustomizeHist(TH1F* hist, Style_t style, Color_t color){
  hist->SetMarkerStyle(style);
  hist-> SetMarkerColor(color);
  hist-> SetLineColor(color);
  hist->SetMarkerSize(1.0);
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
//===================================
// Function to draw histograms in a canva and save it in pdf
// File vector must have signal file first and then all bkg files
// Signal is overlapped to the stack of backgrounds
void DrawHistos(std::vector<TFile*> file,const char *histname,std::vector<Double_t> weights,std::vector<Color_t> color,TString XTitle,TString YTitle,std::array<Double_t,4> leglimits,TString legtitle,std::vector<TString> legentry,const char *outname, TString outfolder){

   // Control size of inputs
   if( !( (file.size() == weights.size()) && (color.size() == legentry.size()) && (file.size() == color.size()) )  ){
    cout << "Error in DrawHistos: number of histograms not defined in " << histname << endl;
    cout << "file size: " << file.size() << "   weights size: " << weights.size() << endl;
    cout << "color size: " << color.size() << "   legend entry size: " << legentry.size() << endl;
    exit(EXIT_FAILURE);
  }

  // Take histograms from files
  std::vector<TH1F*> hist;

  for(Int_t i = 0; i < file.size(); ++i){
    TH1F *h = (TH1F*)file.at(i)->Get(histname)->Clone( (histname + std::to_string(i)).c_str() );
    hist.push_back(h);
  }

  // Scale histograms by Integrated Luminosity, Cross-section and number of generated events
  for(Int_t i = 0; i < hist.size(); ++i){
    hist.at(i)->Sumw2();
    hist.at(i)->Scale(weights.at(i));
    cout << histname << " number " << i << " scaled by " << weights.at(i) << endl;
  }

  // Scale signal to have same area as sum of bkgs
  Double_t BkgIntegral = 0.;

  for(Int_t i = 1; i < hist.size(); ++i){
    BkgIntegral += hist.at(i)-> Integral();
  }

  hist.at(0)->Scale( BkgIntegral/hist.at(0)->Integral() );

  // Write on file
  for(Int_t i = 0; i < hist.size(); ++i){
    hist.at(i)->Write("",TObject::kOverwrite);
  }

  // Canva to draw histograms
  TString canva;
  canva.Form("c%s",histname);

  TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();
  //gPad->SetLogx();

  // Customize signal histo
  CustomizeHist(hist.at(0), kFullCircle, color.at(0));

  // Customize and add bkg histos to stack
  THStack *h_stack = new THStack();
  for(Int_t i = 1; i < hist.size(); ++i){
    hist.at(i)->SetFillColor(color.at(i));
    hist.at(i)-> SetLineColor(kBlue);
    h_stack->Add(hist.at(i));
  }
  
  // Draw histos
  h_stack->Draw("E, HIST");
  hist.at(0)->Draw("E,HIST same");

  // Titles of axes
  h_stack->GetXaxis()->SetTitleSize(.045);
  h_stack->GetYaxis()->SetTitleSize(.045);
  h_stack->GetXaxis()->SetTitle(XTitle);
  h_stack->GetYaxis()->SetTitle(YTitle);

  // Set Range
  Double_t Max = 5e10;
  Double_t Min = 5e0;

  h_stack->SetMaximum(Max);
  h_stack->SetMinimum(Min);
  hist.at(0)->SetMaximum(Max);
  hist.at(0)->SetMinimum(Min);

  // Legend
  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = "#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L=100fb^{-1}}}{" + legtitle + "}";

  legend->SetHeader(head,"C");

  for(Int_t i = 0; i < hist.size(); ++i){
    legend->AddEntry(hist.at(i),legentry.at(i));
  }

  CustomizeLeg(legend);
  legend->Draw();

  TString outpdf;
  outpdf.Form(outfolder + "%s.pdf",outname);
  c->SaveAs(outpdf);

  //c->Close();
}
//===================================
// Function to draw histograms in a canva and save it in pdf
// File vector must have signal file first and then all bkg files
// Signal is stacked with backgrounds
void DrawHistosStack(std::vector<TFile*> file,const char *histname,std::vector<Double_t> weights,std::vector<Color_t> color,TString XTitle,TString YTitle,std::array<Double_t,4> leglimits,TString legtitle,std::vector<TString> legentry,const char *outname, TString outfolder){

   // Control size of inputs
   if( !( (file.size() == weights.size()) && (color.size() == legentry.size()) && (file.size() == color.size()) )  ){
    cout << "Error in DrawHistos: number of histograms not defined in " << histname << endl;
    cout << "file size: " << file.size() << "   weights size: " << weights.size() << endl;
    cout << "color size: " << color.size() << "   legend entry size: " << legentry.size() << endl;
    exit(EXIT_FAILURE);
  }

  // Take histograms from files
  std::vector<TH1F*> hist;

  for(Int_t i = 0; i < file.size(); ++i){
    TH1F *h = (TH1F*)file.at(i)->Get(histname)->Clone( (histname + std::to_string(i)).c_str() );
    hist.push_back(h);
  }

  // Scale histograms by Integrated Luminosity, Cross-section and number of generated events
  for(Int_t i = 0; i < hist.size(); ++i){
    hist.at(i)->Sumw2();
    hist.at(i)->Scale(weights.at(i));
    cout << histname << " number " << i << " scaled by " << weights.at(i) << endl;
  }

  // Scale signal to have same area as sum of bkgs
  Double_t BkgIntegral = 0.;

  for(Int_t i = 1; i < hist.size(); ++i){
    BkgIntegral += hist.at(i)-> Integral();
  }

  hist.at(0)->Scale( BkgIntegral/hist.at(0)->Integral() );

  // Canva to draw histograms
  TString canva;
  canva.Form("cstack%s",histname);

  TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();
  //gPad->SetLogx();

  // Customize and add bkg histos to stack
  THStack *h_stack = new THStack();
  for(Int_t i = 1; i < hist.size(); ++i){
    hist.at(i)->SetFillColor(color.at(i));
    hist.at(i)-> SetLineColor(kBlue);
    h_stack->Add(hist.at(i));
  }

  // Add signal on top of stack
  hist.at(0)->SetFillColor(color.at(0));
  hist.at(0)-> SetLineColor(kBlue);
  h_stack->Add(hist.at(0));
  
  // Draw histos
  h_stack->Draw("E, HIST");

  // Titles of axes
  h_stack->GetXaxis()->SetTitleSize(.045);
  h_stack->GetYaxis()->SetTitleSize(.045);
  h_stack->GetXaxis()->SetTitle(XTitle);
  h_stack->GetYaxis()->SetTitle(YTitle);

  // Set Range
  Double_t Max = 5e10;
  Double_t Min = 5e0;

  h_stack->SetMaximum(Max);
  h_stack->SetMinimum(Min);

  // Legend
  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = "#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L=100fb^{-1}}}{" + legtitle + "}";

  legend->SetHeader(head,"C");

  for(Int_t i = 0; i < hist.size(); ++i){
    legend->AddEntry(hist.at(i),legentry.at(i));
  }

  CustomizeLeg(legend);
  legend->Draw();

  TString outpdf;
  outpdf.Form(outfolder + "%s_stack.pdf",outname);
  c->SaveAs(outpdf);

  //c->Close();
}
//===============================
// Function to draw 2D histograms
void DrawHistos2D(std::vector<TFile*> file,const char *histname,std::vector<Double_t> weights,TString XTitle,TString YTitle,TString ZTitle,std::array<Double_t,4> leglimits,TString legtitle,std::vector<TString> legentry,const char *outname, TString outfolder){

   // Control size of inputs
   if( !( (file.size() == weights.size()) && (file.size() == legentry.size()) )  ){
    cout << "Error in DrawHistos: number of histograms not defined in " << histname << endl;
    cout << "file size: " << file.size() << "   weights size: " << weights.size() << endl;
    cout << "legend entry size: " << legentry.size() << endl;
    exit(EXIT_FAILURE);
  }

  // Take histograms from files
  std::vector<TH2F*> hist;

  for(Int_t i = 0; i < file.size(); ++i){
    TH2F *h = (TH2F*)file.at(i)->Get(histname)->Clone( (histname + std::to_string(i)).c_str() );
    hist.push_back(h);
  }

  // Draw signal and backgrounds separately

  TString canva;
  TString head;
  TString outpdf;

  // Maximum and minimum values of histograms
  Double_t Max = 1e5;
  Double_t Min = 1e-1;

  // Scale histograms by Integrated Luminosity, Cross-section and number of generated events
  for(Int_t i = 0; i < hist.size(); ++i){
    hist.at(i)->Sumw2();
    hist.at(i)->Scale(weights.at(i));
    cout << histname << " number " << i << " scaled by " << weights.at(i) << endl;
  }

  // Scale signal to have same area as sum of bkgs
  Double_t BkgIntegral = 0.;

  for(Int_t i = 1; i < hist.size(); ++i){
    BkgIntegral += hist.at(i)-> Integral();
  }

  hist.at(0)->Scale( BkgIntegral/hist.at(0)->Integral() );

  // Write on file
  for(Int_t i = 0; i < hist.size(); ++i){
    hist.at(i)->Write("",TObject::kOverwrite);
  }

  for(Int_t i = 0; i < hist.size(); ++i){

    // Canva to draw histograms
    canva.Form("c%s_%d",histname,i);

    TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
    gStyle->SetOptStat(0);
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.15);
    //gPad->SetLogy();
    //gPad->SetLogx();
    gPad->SetLogz();

    // Draw histos
    hist.at(i)->Draw("COLZ0");
    hist.at(i)->SetMaximum(Max);
    hist.at(i)->SetMinimum(Min);

    // Titles of axes
    hist.at(i)->GetXaxis()->SetTitleSize(.045);
    hist.at(i)->GetYaxis()->SetTitleSize(.045);
    hist.at(i)->GetZaxis()->SetTitleSize(.045);
    hist.at(i)->GetXaxis()->SetTitle(XTitle);
    hist.at(i)->GetYaxis()->SetTitle(YTitle);
    hist.at(i)->GetZaxis()->SetTitle(ZTitle);

    // Legend
    TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

    head = "#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L=100fb^{-1}}}{" + legtitle + "}";

    legend->SetHeader(head,"C");
    CustomizeLeg(legend);

    // Add the right entry to the legend
    legend->AddEntry(hist.at(i),legentry.at(i));
    legend->Draw();

    outpdf.Form(outfolder + "%s_sample%d.pdf",outname,i);
    c->SaveAs(outpdf);

    //c->Close();
  }

}
//===============================
//---------- END FUNCTIONS -------------------


void MergeHistosArea()
{

  gSystem->Load("libDelphes");

  // Import histograms 
  TString sampledir = "/gfsvol01/atlas/giuliac/plots_and_outputs/20250527_Task21/";

  TString signame = sampledir + "signal_EJ.root";
  TString bkgZqname = sampledir + "bkg_Zjets_q.root";
  TString bkgZgname = sampledir + "bkg_Zjets_g.root";
  TString bkgZZname = sampledir + "bkg_ZZ.root";
  TString bkgHZname = sampledir + "bkg_HZ_SM.root";

  TFile *f_sig = new TFile(signame);
  TFile *f_bkgZq = new TFile(bkgZqname);
  TFile *f_bkgZg = new TFile(bkgZgname);
  TFile *f_bkgZZ = new TFile(bkgZZname);
  TFile *f_bkgHZ = new TFile(bkgHZname);

  // Integrated Luminosity [pb^-1] (arbitrary, no real data here)
  Double_t lumi = 100000;

  // Cross-section [pb] 
  Double_t xsec_bkgZq = 1.250e4;
  Double_t xsec_bkgZg = 7.759e3;
  Double_t xsec_bkgZZ = 7.688e-1;
  Double_t xsec_bkgHZ = 6.455e-2;
  Double_t xsec_sig = xsec_bkgZq + xsec_bkgZg +xsec_bkgZZ + xsec_bkgHZ;

  // Number of MC events
  Double_t ntot_sig = 500000;
  Double_t ntot_bkgZq = 2500000;
  Double_t ntot_bkgZg = 2500000;
  Double_t ntot_bkgZZ = 500000;
  Double_t ntot_bkgHZ = 500000;

  // Prepare info to draw histos and draw them
  std::vector<TFile*> infiles = {f_sig,f_bkgHZ,f_bkgZZ,f_bkgZg,f_bkgZq};
  std::vector<Color_t> colors = {kBlack,kAzure+5,kGreen+1,kMagenta-10,kOrange-8};
  //std::vector<TString> legentry = {TString::Format("f #bar{f} -> H Z -> qv1 #bar{qv1} l^{+} l^{-}, #sigma = %.0f pb",xsec_sig),TString::Format("HZ -> b #bar{b} l^{+} l^{-}, #sigma = %.4f pb",xsec_bkgHZ),TString::Format("ZZ, #sigma = %.3f pb",xsec_bkgZZ),TString::Format("Z + jets (g), #sigma = %.0f pb",xsec_bkgZg),TString::Format("Z + jets (q), #sigma = %.0f pb",xsec_bkgZq)};
  std::vector<TString> legentry = {TString::Format("f #bar{f} -> H Z -> qv1 #bar{qv1} l^{+} l^{-}, Area = #Sigma Bkgs"),TString::Format("HZ -> b #bar{b} l^{+} l^{-}, #sigma = %.4f pb",xsec_bkgHZ),TString::Format("ZZ, #sigma = %.3f pb",xsec_bkgZZ),TString::Format("Z + jets (g), #sigma = %.0f pb",xsec_bkgZg),TString::Format("Z + jets (q), #sigma = %.0f pb",xsec_bkgZq)};

  // Calculate weights
  std::vector<Double_t> xsec = {xsec_sig,xsec_bkgHZ,xsec_bkgZZ,xsec_bkgZg,xsec_bkgZq};
  std::vector<Double_t> ntot = {ntot_sig,ntot_bkgHZ,ntot_bkgZZ,ntot_bkgZg,ntot_bkgZq};
  std::vector<Double_t> weights;

  for(Int_t i = 0; i < infiles.size(); ++i){
    weights.push_back( lumi*xsec.at(i) / ntot.at(i) );
  }

  // outnames with "area"
  sampledir += "area_" ;

  // write all scaled histograms on file
  TString outname = sampledir + "MergedHistos.root";

  TFile *f = new TFile(outname,"RECREATE");

  // 1D histograms

  std::array<Double_t,4> leglimitsInvmassS = {0.452381,0.543403,0.763158,0.840278};
  DrawHistos(infiles,"Invariant_mass_of_leading_and_subleading_small_jets",weights,colors,"Invariant Mass [GeV/c^{2}]","#events",leglimitsInvmassS,"Invariant Mass of first 2 Small Jets",legentry,"InvMass_small",sampledir);
  DrawHistosStack(infiles,"Invariant_mass_of_leading_and_subleading_small_jets",weights,colors,"Invariant Mass [GeV/c^{2}]","#events",leglimitsInvmassS,"Invariant Mass of first 2 Small Jets",legentry,"InvMass_small",sampledir);

  std::array<Double_t,4> leglimitsInvmassSA = {0.472431,0.546875,0.738095,0.84375};
  DrawHistos(infiles,"Invariant_mass_of_leading_and_subleading_small_jets_thesis_abstract",weights,colors,"Invariant Mass [GeV/c^{2}]","#events",leglimitsInvmassSA,"Invariant Mass of first 2 Small Jets",legentry,"InvMass_small_abstract",sampledir);

  std::array<Double_t,4> leglimitsInvmassL = {0.452381,0.543403,0.763158,0.840278};
  DrawHistos(infiles,"Invariant_mass_of_leading_and_subleading_large_jets",weights,colors,"Invariant Mass [GeV/c^{2}]","#events",leglimitsInvmassL,"Invariant Mass of first 2 Large Jets",legentry,"InvMass_large",sampledir);
  DrawHistosStack(infiles,"Invariant_mass_of_leading_and_subleading_large_jets",weights,colors,"Invariant Mass [GeV/c^{2}]","#events",leglimitsInvmassL,"Invariant Mass of first 2 Large Jets",legentry,"InvMass_large",sampledir);

  std::array<Double_t,4> leglimitsDeltaR = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"DeltaR_between_leading_and_subleading_small_jets",weights,colors,"#Delta R","# events",leglimitsDeltaR,"#Delta R between first 2 Small Jets",legentry,"DeltaRS",sampledir);
  DrawHistosStack(infiles,"DeltaR_between_leading_and_subleading_small_jets",weights,colors,"#Delta R","# events",leglimitsDeltaR,"#Delta R between first 2 Small Jets",legentry,"DeltaRS",sampledir);

  std::array<Double_t,4> leglimitsDeltaPhi = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"DeltaPhi_between_leading_and_subleading_small_jets",weights,colors,"|#Delta #varphi|","# events",leglimitsDeltaPhi,"|#Delta #varphi| between first 2 Small Jets",legentry,"DeltaPhiS",sampledir);
  DrawHistosStack(infiles,"DeltaPhi_between_leading_and_subleading_small_jets",weights,colors,"|#Delta #varphi|","# events",leglimitsDeltaPhi,"|#Delta #varphi| between first 2 Small Jets",legentry,"DeltaPhiS",sampledir);

  DrawHistos(infiles,"DeltaR_between_leading_and_subleading_large_jets",weights,colors,"#Delta R","# events",leglimitsDeltaR,"#Delta R between first 2 Large Jets",legentry,"DeltaRL",sampledir);
  DrawHistosStack(infiles,"DeltaR_between_leading_and_subleading_large_jets",weights,colors,"#Delta R","# events",leglimitsDeltaR,"#Delta R between first 2 Large Jets",legentry,"DeltaRL",sampledir);

  DrawHistos(infiles,"DeltaPhi_between_leading_and_subleading_large_jets",weights,colors,"|#Delta #varphi|","# events",leglimitsDeltaPhi,"|#Delta #varphi| between first 2 Large Jets",legentry,"DeltaPhiL",sampledir);
  DrawHistosStack(infiles,"DeltaPhi_between_leading_and_subleading_large_jets",weights,colors,"|#Delta #varphi|","# events",leglimitsDeltaPhi,"|#Delta #varphi| between first 2 Large Jets",legentry,"DeltaPhiL",sampledir);

  std::array<Double_t,4> leglimitsPT1S = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"leading_small_jet_pt",weights,colors,"p_{T} [GeV/c]","# events",leglimitsPT1S,"pT of Leading Small Jet",legentry,"JetPT1S",sampledir);
  DrawHistosStack(infiles,"leading_small_jet_pt",weights,colors,"p_{T} [GeV/c]","# events",leglimitsPT1S,"pT of Leading Small Jet",legentry,"JetPT1S",sampledir);

  std::array<Double_t,4> leglimitsPT2S = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"subleading_small_jet_pt",weights,colors,"p_{T} [GeV/c]","# events",leglimitsPT2S,"pT of Subleading Small Jet",legentry,"JetPT2S",sampledir);
  DrawHistosStack(infiles,"subleading_small_jet_pt",weights,colors,"p_{T} [GeV/c]","# events",leglimitsPT2S,"pT of Subleading Small Jet",legentry,"JetPT2S",sampledir);

  std::array<Double_t,4> leglimitsPT1L = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"leading_large_jet_pt",weights,colors,"p_{T} [GeV/c]","# events",leglimitsPT1L,"pT of Leading Large Jet",legentry,"JetPT1L",sampledir);
  DrawHistosStack(infiles,"leading_large_jet_pt",weights,colors,"p_{T} [GeV/c]","# events",leglimitsPT1L,"pT of Leading Large Jet",legentry,"JetPT1L",sampledir);

  std::array<Double_t,4> leglimitsPT2L = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"subleading_large_jet_pt",weights,colors,"p_{T} [GeV/c]","# events",leglimitsPT2L,"pT of Subleading Large Jet",legentry,"JetPT2L",sampledir);
  DrawHistosStack(infiles,"subleading_large_jet_pt",weights,colors,"p_{T} [GeV/c]","# events",leglimitsPT2L,"pT of Subleading Large Jet",legentry,"JetPT2L",sampledir);

  std::array<Double_t,4> leglimitsNC1 = {0.452381,0.543403,0.763158,0.840278};
  DrawHistos(infiles,"number_of_constituents_of_leading_large_jet",weights,colors,"Number of Constituents","# events",leglimitsNC1,"Number of Constituents of Leading Large Jet",legentry,"NConst1",sampledir);
  DrawHistosStack(infiles,"number_of_constituents_of_leading_large_jet",weights,colors,"Number of Constituents","# events",leglimitsNC1,"Number of Constituents of Leading Large Jet",legentry,"NConst1",sampledir);

  std::array<Double_t,4> leglimitsNC2 = {0.452381,0.543403,0.763158,0.840278};
  DrawHistos(infiles,"number_of_constituents_of_subleading_large_jet",weights,colors,"Number of Constituents","# events",leglimitsNC2,"Number of Constituents of Subleading Large Jet",legentry,"NConst2",sampledir);
  DrawHistosStack(infiles,"number_of_constituents_of_subleading_large_jet",weights,colors,"Number of Constituents","# events",leglimitsNC2,"Number of Constituents of Subleading Large Jet",legentry,"NConst2",sampledir);

  std::array<Double_t,4> leglimitsNS = {0.452381,0.543403,0.763158,0.840278};
  DrawHistos(infiles,"number_of_small_jets",weights,colors,"Number of Jets","# events",leglimitsNS,"Number of Small Jets",legentry,"Nsmall",sampledir);
  DrawHistosStack(infiles,"number_of_small_jets",weights,colors,"Number of Jets","# events",leglimitsNS,"Number of Small Jets",legentry,"Nsmall",sampledir);

  std::array<Double_t,4> leglimitsNbS = {0.452381,0.543403,0.763158,0.840278};
  DrawHistos(infiles,"number_of_btagged_small_jets",weights,colors,"Number of Jets","# events",leglimitsNbS,"Number of b-tagged Small Jets",legentry,"Nbsmall",sampledir);
  DrawHistosStack(infiles,"number_of_btagged_small_jets",weights,colors,"Number of Jets","# events",leglimitsNbS,"Number of b-tagged Small Jets",legentry,"Nbsmall",sampledir);

  std::array<Double_t,4> leglimitsNL = {0.452381,0.543403,0.763158,0.840278};
  DrawHistos(infiles,"number_of_large_jets",weights,colors,"Number of Jets","# events",leglimitsNL,"Number of Large Jets",legentry,"Nlarge",sampledir);
  DrawHistosStack(infiles,"number_of_large_jets",weights,colors,"Number of Jets","# events",leglimitsNL,"Number of Large Jets",legentry,"Nlarge",sampledir);

  std::array<Double_t,4> leglimitsPTF1L = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"Prompt_Track_Fraction_of_leading_large_jet",weights,colors,"PTF","# events",leglimitsPTF1L,"PTF of Leading Large Jet",legentry,"PTF1L",sampledir);
  DrawHistosStack(infiles,"Prompt_Track_Fraction_of_leading_large_jet",weights,colors,"PTF","# events",leglimitsPTF1L,"PTF of Leading Large Jet",legentry,"PTF1L",sampledir);

  std::array<Double_t,4> leglimitsPTF2L = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"Prompt_Track_Fraction_of_subleading_large_jet",weights,colors,"PTF","# events",leglimitsPTF2L,"PTF of Subleading Large Jet",legentry,"PTF2L",sampledir);
  DrawHistosStack(infiles,"Prompt_Track_Fraction_of_subleading_large_jet",weights,colors,"PTF","# events",leglimitsPTF2L,"PTF of Subleading Large Jet",legentry,"PTF2L",sampledir);

  std::array<Double_t,4> leglimitsPTF1S = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"Prompt_Track_Fraction_of_leading_small_jet",weights,colors,"PTF","# events",leglimitsPTF1S,"PTF of Leading Small Jet",legentry,"PTF1S",sampledir);
  DrawHistosStack(infiles,"Prompt_Track_Fraction_of_leading_small_jet",weights,colors,"PTF","# events",leglimitsPTF1S,"PTF of Leading Small Jet",legentry,"PTF1S",sampledir);

  std::array<Double_t,4> leglimitsPTF2S = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"Prompt_Track_Fraction_of_subleading_small_jet",weights,colors,"PTF","# events",leglimitsPTF2S,"PTF of Subleading Small Jet",legentry,"PTF2S",sampledir);
  DrawHistosStack(infiles,"Prompt_Track_Fraction_of_subleading_small_jet",weights,colors,"PTF","# events",leglimitsPTF2S,"PTF of Subleading Small Jet",legentry,"PTF2S",sampledir);

  std::array<Double_t,4> leglimitsPTF1bS = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"Prompt_Track_Fraction_of_leading_btagged_small_jet",weights,colors,"PTF","# events",leglimitsPTF1bS,"PTF of Leading b-tagged Small Jet",legentry,"PTF1bS",sampledir);
  DrawHistosStack(infiles,"Prompt_Track_Fraction_of_leading_btagged_small_jet",weights,colors,"PTF","# events",leglimitsPTF1bS,"PTF of Leading b-tagged Small Jet",legentry,"PTF1bS",sampledir);

  std::array<Double_t,4> leglimitsECF1L = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"2point_Energy_Correlation_div_pT_of_leading_large_jet",weights,colors,"ECF2/pT [GeV/c]","# events",leglimitsECF1L,"ECF2/pT of Leading Large Jet",legentry,"ECF1L",sampledir);
  DrawHistosStack(infiles,"2point_Energy_Correlation_div_pT_of_leading_large_jet",weights,colors,"ECF2/pT [GeV/c]","# events",leglimitsECF1L,"ECF2/pT of Leading Large Jet",legentry,"ECF1L",sampledir);

  std::array<Double_t,4> leglimitsECF2L = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"2point_Energy_Correlation_div_pT_of_subleading_large_jet",weights,colors,"ECF2/pT [GeV/c]","# events",leglimitsECF2L,"ECF2/pT of Subleading Large Jet",legentry,"ECF2L",sampledir);
  DrawHistosStack(infiles,"2point_Energy_Correlation_div_pT_of_subleading_large_jet",weights,colors,"ECF2/pT [GeV/c]","# events",leglimitsECF2L,"ECF2/pT of Subleading Large Jet",legentry,"ECF2L",sampledir);

  std::array<Double_t,4> leglimitsLepPT = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"pt_of_first_two_leptons",weights,colors,"p_{T} [GeV/c]","# events",leglimitsLepPT,"pT of Z tagging Leptons",legentry,"Lep2PT",sampledir);
  DrawHistosStack(infiles,"pt_of_first_two_leptons",weights,colors,"p_{T} [GeV/c]","# events",leglimitsLepPT,"pT of Z tagging Leptons",legentry,"Lep2PT",sampledir);

  std::array<Double_t,4> leglimitsLepR = {0.452381,0.571181,0.763158,0.838542};
  DrawHistos(infiles,"DeltaR_between_first_two_leptons",weights,colors,"#Delta R","# events",leglimitsLepR,"#Delta R between Z tagging Leptons",legentry,"LepDeltaR",sampledir);
  DrawHistosStack(infiles,"DeltaR_between_first_two_leptons",weights,colors,"#Delta R","# events",leglimitsLepR,"#Delta R between Z tagging Leptons",legentry,"LepDeltaR",sampledir);

  // 2D histograms

  std::array<Double_t,4> leglimitsDRP = {0.362155,0.65625,0.672932,0.878472};
  DrawHistos2D(infiles,"DeltaRsmall_vs_DeltaPhismall",weights,"#Delta R","|#Delta #varphi|","# events",leglimitsDRP,"#Delta R vs |#Delta #varphi| between first 2 Small Jets",legentry,"DeltaRs_DeltaPhis",sampledir);

  std::array<Double_t,4> leglimitsDRR = {0.362155,0.65625,0.672932,0.878472};
  DrawHistos2D(infiles,"DeltaRsmall_vs_DeltaRlarge",weights,"#Delta R between first 2 Small Jets","#Delta R between first 2 Large Jets","# events",leglimitsDRR,"#Delta R between first 2 Small Jets vs #Delta R between first 2 Large Jets",legentry,"DeltaRs_DeltaRl",sampledir);

  std::array<Double_t,4> leglimitsDPP = {0.362155,0.65625,0.672932,0.878472};
  DrawHistos2D(infiles,"DeltaPhismall_vs_DeltaPhilarge",weights,"|#Delta #varphi| between first 2 Small Jets","|#Delta #varphi| between first 2 Large Jets","# events",leglimitsDPP,"|#Delta #varphi| between first 2 Small Jets vs |#Delta #varphi| between first 2 Large Jets",legentry,"DeltaPhis_DeltaPhil",sampledir);

  std::array<Double_t,4> leglimitsPT1s_PT2s = {0.235589,0.689236,0.546366,0.869792};
  DrawHistos2D(infiles,"JetPT1small_vs_JetPT2small",weights,"p_{T} of Leading Small Jet [GeV/c]","p_{T} of Subleading Small Jet [GeV/c]","# events",leglimitsPT1s_PT2s,"L Small Jet pT vs SubL Small Jet pT",legentry,"JetPT1s_JetPT2s",sampledir);

  std::array<Double_t,4> leglimitsPT1l_PT2l = {0.235589,0.689236,0.546366,0.869792};
  DrawHistos2D(infiles,"JetPT1large_vs_JetPT2large",weights,"p_{T} of Leading Large Jet [GeV/c]","p_{T} of Subleading Large Jet [GeV/c]","# events",leglimitsPT1l_PT2l,"L Large Jet pT vs SubL Large Jet pT",legentry,"JetPT1l_JetPT2l",sampledir);

  std::array<Double_t,4> leglimitsIMs_IMl = {0.408521,0.119792,0.697995,0.28125};
  DrawHistos2D(infiles,"IMsmall_vs_IMlarge",weights,"Invariant Mass of first 2 Small Jets [GeV/c^2]","Invariant Mass of first 2 Large Jets [GeV/c^2]","# events",leglimitsIMs_IMl,"IM of first 2 Small Jets vs IM of first 2 Large Jets",legentry,"IMs_IMl",sampledir);

  std::array<Double_t,4> leglimitsIMs_DRs = {0.235589,0.689236,0.546366,0.869792};
  DrawHistos2D(infiles,"IMsmall_vs_DeltaRsmall",weights,"Invariant Mass of first 2 Small Jets [GeV/c^2]","#Delta R between first 2 Small Jets","# events",leglimitsIMs_DRs,"IM vs #Delta R of first 2 Small Jets",legentry,"IMs_DRs",sampledir);

  std::array<Double_t,4> leglimitsIMl_DRl = {0.235589,0.689236,0.546366,0.869792};
  DrawHistos2D(infiles,"IMlarge_vs_DeltaRlarge",weights,"Invariant Mass of first 2 Large Jets [GeV/c^2]","#Delta R between first 2 Large Jets","# events",leglimitsIMl_DRl,"IM vs #Delta R of first 2 Large Jets",legentry,"IMl_DRl",sampledir);

}//end MergeHistos
