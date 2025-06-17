/*
Author: Giulia Cossutti

Macro to merge histograms for signal and backgrounds.
Bkgs are scaled with pythia's calculated cross-sections or measured cross-sections, 
Signal is scaled to have same area as sum of bkgs.
Histograms are saved in files.
Events are divided in cathegories: A, B.

From inside the /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes directory run with
root -l /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250617_Task23/MergeHistosArea.C 

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
  Double_t Max = 5e8;
  Double_t Min = 5e0;

  h_stack->SetMaximum(Max);
  h_stack->SetMinimum(Min);
  hist.at(0)->SetMaximum(Max);
  hist.at(0)->SetMinimum(Min);

  // Legend
  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, Cathegory %.1s}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L=100fb^{-1}}}{" + legtitle + "}", histname);

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
  Double_t Max = 5e8;
  Double_t Min = 5e0;

  h_stack->SetMaximum(Max);
  h_stack->SetMinimum(Min);

  // Legend
  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, Cathegory %.1s}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L=100fb^{-1}}}{" + legtitle + "}", histname);

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
  Double_t Max = 1e6;
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

    head = TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, Cathegory %.1s}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L=100fb^{-1}}}{" + legtitle + "}", histname);

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
  TString sampledir = "/gfsvol01/atlas/giuliac/plots_and_outputs/20250617_Task23/";

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
  //Double_t xsec_sig = xsec_bkgZq + xsec_bkgZg;

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

  // Event cathegories
  std::vector<TString> cat = {"A_","B_"};

  TFile *f = new TFile(outname,"RECREATE");

  for(Int_t i = 0; i < cat.size(); ++i){

    // 1D histograms

    std::array<Double_t,4> leglimitsDeltaR = {0.472431,0.543403,0.783208,0.817708};
    DrawHistos(infiles,cat.at(i) + "DeltaR_between_leading_and_subleading_small_jets",weights,colors,"#Delta R","# events",leglimitsDeltaR,"#Delta R between first 2 Small Jets",legentry,cat.at(i) + "DeltaRS",sampledir);

    std::array<Double_t,4> leglimitsNbS = {0.472431,0.543403,0.783208,0.817708};
    DrawHistos(infiles,cat.at(i) + "number_of_btagged_small_jets",weights,colors,"Number of Jets","# events",leglimitsNbS,"Number of b-tagged Small Jets",legentry,cat.at(i) + "Nbsmall",sampledir);

    std::array<Double_t,4> leglimitsLepPT = {0.472431,0.543403,0.783208,0.817708};
    DrawHistos(infiles,cat.at(i) + "pt_of_first_two_leptons",weights,colors,"p_{T} [GeV/c]","# events",leglimitsLepPT,"pT of Z tagging Leptons",legentry,cat.at(i) + "Lep2PT",sampledir);

    // 2D histograms

    std::array<Double_t,4> leglimitsNs_Nbs = {0.235589,0.689236,0.546366,0.869792};
    DrawHistos2D(infiles,cat.at(i) + "Nsmall_vs_Nbsmall",weights,"Number of Small Jets","Number of b-tagged Small Jets","# events",leglimitsNs_Nbs,"Number of Small vs b-tagged Small Jets",legentry,cat.at(i) + "Ns_Nbs",sampledir);

  }

}//end MergeHistos
