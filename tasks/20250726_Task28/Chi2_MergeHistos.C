/*
Author: Giulia Cossutti

Macro to merge histograms for signal and backgrounds.
Bkgs are scaled with pythia's calculated cross-sections or measured cross-sections, 
Signal is scaled to have same area as sum of bkgs.
Histograms are saved in files.

Calculation of Chi2 with different signal xsections to determine 95% CL limits.
Both statistical and systematic uncertainties are considered.
Chi2 graphs with single systematics, plots of single systematics.

Events are divided in categories: A, B.
Pre-selection is applied.
Cut1 is applied (number of b-tagged small jets >= 1) for both A and B events.
Cut2 is applied:
- A: DeltaR between first 2 small jets < 2.2
- B: DeltaR between first 2 small jets < 4.8
Cut3 is applied:
- A: pT of sum of 2 leptons by Z decay > 96 GeV/c
Cut4 is applied:
- A: DeltaR between first 2 small jets <= -0.09 * mjj of first 2 small jets + 8.5
- B: ECF2/pT of Leading Large Jet > 2 GeV/c

From inside the /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes directory run with
root -l /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250726_Task28/Chi2_MergeHistos.C 

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
// Function to customize Chi2vsxsec_sig graphs:
void CustomizeGr(TGraph* cl_gr, Color_t color){
  cl_gr->SetMarkerStyle(kFullCircle);
  cl_gr->SetMarkerColor(kBlack);
  cl_gr->SetLineColor(color);
  cl_gr->SetLineWidth(2);
  cl_gr->GetXaxis()->SetTitleSize(.045);
  cl_gr->GetYaxis()->SetTitleSize(.045);
  cl_gr->GetXaxis()->SetTitle("#sigma_{sig} [pb]");
  cl_gr->GetYaxis()->SetTitle("#chi^{2}");
}
//===================================
// Function to draw histograms in a canva and save it in pdf
// File vector must have signal file first and then all bkg files
// Signal is overlapped to the stack of backgrounds
void DrawHistos(Bool_t area,std::vector<TFile*> file,const char *histname,std::vector<Double_t> weights,std::vector<Color_t> color,TString XTitle,TString YTitle,std::array<Double_t,4> leglimits,TString legtitle,std::vector<TString> legentry,const char *outname, TString outfolder){

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
  // Set Bin Error = sqrt(Bin Content)
  for(Int_t i = 0; i < hist.size(); ++i){
    hist.at(i)->Sumw2();
    hist.at(i)->Scale(weights.at(i));
    cout << histname << " number " << i << " scaled by " << weights.at(i) << endl;

    for(Int_t k = 1; k <= hist.at(i)->GetXaxis()->GetNbins(); ++k){
      hist.at(i)->SetBinError(k, sqrt(hist.at(i)->GetBinContent(k)) );
    }
  }

  if(area == kTRUE){
    // Scale signal to have same area as sum of bkgs
    Double_t BkgIntegral = 0.;

    for(Int_t i = 1; i < hist.size(); ++i){
      BkgIntegral += hist.at(i)-> Integral();
    }

    hist.at(0)->Scale( BkgIntegral/hist.at(0)->Integral() );
    for(Int_t k = 1; k <= hist.at(0)->GetXaxis()->GetNbins(); ++k){
      hist.at(0)->SetBinError(k, sqrt(hist.at(0)->GetBinContent(k)) );
    }
  }

  // Write on file
  for(Int_t i = 0; i < hist.size(); ++i){
    hist.at(i)->Write("",TObject::kSingleKey);
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
  Double_t Max = 5e7;
  Double_t Min = 5e-1;

  h_stack->SetMaximum(Max);
  h_stack->SetMinimum(Min);
  hist.at(0)->SetMaximum(Max);
  hist.at(0)->SetMinimum(Min);

  // Legend
  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, Category %.1s}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L=360fb^{-1}}}{" + legtitle + "}", outname);

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
void DrawHistosStack(Bool_t area, std::vector<TFile*> file,const char *histname,std::vector<Double_t> weights,std::vector<Color_t> color,TString XTitle,TString YTitle,std::array<Double_t,4> leglimits,TString legtitle,std::vector<TString> legentry,const char *outname, TString outfolder){

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
  // Set Bin Error = sqrt(Bin Content)
  for(Int_t i = 0; i < hist.size(); ++i){
    hist.at(i)->Sumw2();
    hist.at(i)->Scale(weights.at(i));
    cout << histname << " number " << i << " scaled by " << weights.at(i) << endl;

    for(Int_t k = 1; k <= hist.at(i)->GetXaxis()->GetNbins(); ++k){
      hist.at(i)->SetBinError(k, sqrt(hist.at(i)->GetBinContent(k)) );
    }
  }

  if(area == kTRUE){
    // Scale signal to have same area as sum of bkgs
    Double_t BkgIntegral = 0.;

    for(Int_t i = 1; i < hist.size(); ++i){
      BkgIntegral += hist.at(i)-> Integral();
    }

    hist.at(0)->Scale( BkgIntegral/hist.at(0)->Integral() );
    for(Int_t k = 1; k <= hist.at(0)->GetXaxis()->GetNbins(); ++k){
      hist.at(0)->SetBinError(k, sqrt(hist.at(0)->GetBinContent(k)) );
    }
  }

  // Write on file
  for(Int_t i = 0; i < hist.size(); ++i){
    hist.at(i)->Write("",TObject::kSingleKey);
  }

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
  Double_t Max = 5e7;
  Double_t Min = 5e-1;

  h_stack->SetMaximum(Max);
  h_stack->SetMinimum(Min);

  // Legend
  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, Category %.1s}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L=360fb^{-1}}}{" + legtitle + "}", outname);

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
void DrawHistos2D(Bool_t area,std::vector<TFile*> file,const char *histname,std::vector<Double_t> weights,TString XTitle,TString YTitle,TString ZTitle,std::array<Double_t,4> leglimits,TString legtitle,std::vector<TString> legentry,const char *outname, TString outfolder){

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
  Double_t Max = 1e4;
  Double_t Min = 1e-1;

  // Scale histograms by Integrated Luminosity, Cross-section and number of generated events
  // Set Bin Error = sqrt(Bin Content)
  for(Int_t i = 0; i < hist.size(); ++i){
    hist.at(i)->Sumw2();
    hist.at(i)->Scale(weights.at(i));
    cout << histname << " number " << i << " scaled by " << weights.at(i) << endl;

    for(Int_t kx = 1; kx <= hist.at(i)->GetXaxis()->GetNbins(); ++kx){
      for(Int_t ky = 1; ky <= hist.at(i)->GetYaxis()->GetNbins(); ++ky){
        hist.at(i)->SetBinError(kx,ky, sqrt(hist.at(i)->GetBinContent(kx,ky)) );
      }
    }

  }

  if(area == kTRUE){
    // Scale signal to have same area as sum of bkgs
    Double_t BkgIntegral = 0.;

    for(Int_t i = 1; i < hist.size(); ++i){
      BkgIntegral += hist.at(i)-> Integral();
    }

    hist.at(0)->Scale( BkgIntegral/hist.at(0)->Integral() );
    for(Int_t kx = 1; kx <= hist.at(0)->GetXaxis()->GetNbins(); ++kx){
      for(Int_t ky = 1; ky <= hist.at(0)->GetYaxis()->GetNbins(); ++ky){
        hist.at(0)->SetBinError(kx,ky, sqrt(hist.at(0)->GetBinContent(kx,ky)) );
      }
    }
  }

  // Write on file
  for(Int_t i = 0; i < hist.size(); ++i){
    hist.at(i)->Write("",TObject::kSingleKey);
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

    head = TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, Category %.1s}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L=360fb^{-1}}}{" + legtitle + "}", outname);

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
// Function to draw 2D histograms of 1sigma systematics
void DrawSys(TFile* file,const char *histname,TString XTitle,TString YTitle,TString ZTitle,std::array<Double_t,4> leglimits,TString legtitle,TString legentry,const char *outname, TString outfolder){

  // Take histograms from files
  TH2F *hist = (TH2F*)file->Get(histname)->Clone( histname );

  TString canva;
  TString head;
  TString outpdf;

  // Maximum and minimum values of histograms
  Double_t Max = 1e4;
  Double_t Min = 1e-1;

  // Canva to draw histograms
  canva.Form("csys%s",histname);

  TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  //gPad->SetLogy();
  //gPad->SetLogx();
  gPad->SetLogz();

  // Draw histos
  hist->Draw("COLZ0");
  hist->SetMaximum(Max);
  hist->SetMinimum(Min);

  // Titles of axes
  hist->GetXaxis()->SetTitleSize(.045);
  hist->GetYaxis()->SetTitleSize(.045);
  hist->GetZaxis()->SetTitleSize(.045);
  hist->GetXaxis()->SetTitle(XTitle);
  hist->GetYaxis()->SetTitle(YTitle);
  hist->GetZaxis()->SetTitle(ZTitle);

  // Legend
  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  head = TString::Format("#splitline{Category %.1s, PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L=360fb^{-1}}{" + legtitle + "}", outname);

  legend->SetHeader(head,"C");
  CustomizeLeg(legend);

  // Add the right entry to the legend
  legend->AddEntry(hist,legentry);
  legend->Draw();

  outpdf.Form(outfolder + "%s.pdf",outname);
  c->SaveAs(outpdf);

  // Write on file
  hist->Write("",TObject::kSingleKey);

  //c->Close();

}
//===================================
// Function to calculate and plot Chi2 and CL limits
// File vector must have signal file first and then all bkg files
void Chi2CL(std::vector<TFile*> file,const char *histname,std::vector<Double_t> xsec_sig,std::vector<Double_t> xsec,std::vector<Double_t> ntot,Double_t lumi,std::array<Double_t,4> leglimits,TString legtitle,const char *outname, TString outfolder){

   // Control size of inputs
   if( !( (file.size() == xsec.size()) && (file.size() == ntot.size()) )  ){
    cout << "Error in Chi2CL: number of histograms not defined in " << histname << endl;
    cout << "file size: " << file.size() << "   xsec size: " << xsec.size() << endl;
    cout << "ntot size: " << ntot.size() << endl;
    exit(EXIT_FAILURE);
  }

  // Take histograms from files
  std::vector<TH1F*> hist;

  for(Int_t i = 0; i < file.size(); ++i){
    TH1F *h = (TH1F*)file.at(i)->Get(histname)->Clone( (histname + std::to_string(i)).c_str() );
    hist.push_back(h);
  }

  // Scale histograms by Integrated Luminosity, Cross-section and number of generated events
  std::vector<Double_t> weights;

  for(Int_t i = 0; i < file.size(); ++i){
    weights.push_back( lumi*xsec.at(i) / ntot.at(i) );
  }

  for(Int_t i = 1; i < hist.size(); ++i){
    hist.at(i)->Sumw2();
    hist.at(i)->Scale(weights.at(i));
  }

  // Calculate Chi2 for each value of xsec_sig

  // Vector of signal histograms with different xsec_sig
  std::vector<TH1F*> hist_sig;

  for(Int_t j = 0; j < xsec_sig.size(); ++j){
    TH1F *h = (TH1F*)file.at(0)->Get(histname)->Clone( (histname + std::to_string( j + file.size()) ).c_str() );
    hist_sig.push_back(h);
  }

  // Vector of Chi2 with different xsec_sig
  Double_t Chi2j = 0.;
  std::vector<Double_t> Chi2;

  // Histogram of sum of bkgs
  TH1F *hist_bkg = (TH1F*)hist.at(1)->Clone( (histname +  std::to_string( xsec_sig.size() + file.size()) ).c_str() );
  for(Int_t i = 2; i < hist.size(); ++i){
    hist_bkg->Add(hist.at(i));
  }

  for(Int_t j = 0; j < xsec_sig.size(); ++j){

    // Scale properly the signal
    weights.at(0) = lumi*xsec_sig.at(j) / ntot.at(0);
    hist_sig.at(j)->Sumw2();
    hist_sig.at(j)->Scale(weights.at(0));

    // Calculate Chi2
    Chi2j = 0.;

    for(Int_t k = 1; k <= hist_sig.at(j)->GetXaxis()->GetNbins(); ++k){
      Chi2j += pow( hist_sig.at(j)->GetBinContent(k) ,2) / hist_bkg->GetBinContent(k) ;
    }

    Chi2.push_back(Chi2j);
  }

  // Graph of Chi2 vs xsec_sig

  TGraph *cl_gr = new TGraph();
  cl_gr->SetName(TString::Format("CLgraph_%s", histname));
  cl_gr->SetMarkerStyle(kFullCircle);
  cl_gr->SetMarkerColor(kBlack);
  cl_gr->SetLineColor(kGreen+2);
  cl_gr->SetLineWidth(2);
  cl_gr->GetXaxis()->SetTitleSize(.045);
  cl_gr->GetYaxis()->SetTitleSize(.045);
  cl_gr->GetXaxis()->SetTitle("#sigma_{sig} [pb]");
  cl_gr->GetYaxis()->SetTitle("#chi^{2}");

  for(Int_t j = 0; j < xsec_sig.size(); ++j){
    cl_gr->AddPoint(xsec_sig.at(j),Chi2.at(j));
  }

  // CL limit

  Double_t CL = TMath::ChisquareQuantile(0.95, hist_bkg->GetXaxis()->GetNbins());

  TLine* hor = new TLine(xsec_sig.at(0) , CL, xsec_sig.at(xsec_sig.size()-1), CL);
  hor->SetLineColor(kRed);
  hor->SetLineWidth(2);

  TString canva;
  canva.Form("cCL%s",histname);

  TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  //gPad->SetLogy();
 
  cl_gr->Draw("APL");
  hor->Draw("L same");

  // Legend
  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, Category %.1s}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L= %.0f fb^{-1}}}{ 95%% CL limits on #sigma_{sig} using " + legtitle + "}", outname, lumi/1000.);

  legend->SetHeader(head,"C");
  CustomizeLeg(legend);

  // Add the right entry to the legend
  legend->AddEntry(hor,"95% CL");
  legend->Draw();
 
  TString outpdf;
  outpdf.Form(outfolder + "CL_%s.pdf",outname);
  c->SaveAs(outpdf);

  //c->Close();

  // Write graph on file
  cl_gr->Write("",TObject::kSingleKey);
}
//===================================
// Function to calculate and plot Chi2 and CL limits
// File vector must have signal file first and then all bkg files
void Chi2CL2D(std::vector<TFile*> file,const char *histname,std::vector<Double_t> xsec_sig,std::vector<Double_t> xsec,std::vector<Double_t> ntot,Double_t lumi,std::array<Double_t,4> leglimits,TString legtitle,const char *outname, TString outfolder){

   // Control size of inputs
   if( !( (file.size() == xsec.size()) && (file.size() == ntot.size()) )  ){
    cout << "Error in Chi2CL2D: number of histograms not defined in " << histname << endl;
    cout << "file size: " << file.size() << "   xsec size: " << xsec.size() << endl;
    cout << "ntot size: " << ntot.size() << endl;
    exit(EXIT_FAILURE);
  }

  // Take histograms from files
  std::vector<TH2F*> hist;

  for(Int_t i = 0; i < file.size(); ++i){
    TH2F *h = (TH2F*)file.at(i)->Get(histname)->Clone( (histname + std::to_string(i)).c_str() );
    hist.push_back(h);
  }

  // Scale histograms by Integrated Luminosity, Cross-section and number of generated events
  std::vector<Double_t> weights;

  for(Int_t i = 0; i < file.size(); ++i){
    weights.push_back( lumi*xsec.at(i) / ntot.at(i) );
  }

  for(Int_t i = 1; i < hist.size(); ++i){
    hist.at(i)->Sumw2();
    hist.at(i)->Scale(weights.at(i));
  }

  // Calculate Chi2 for each value of xsec_sig

  // Vector of signal histograms with different xsec_sig
  std::vector<TH2F*> hist_sig;

  for(Int_t j = 0; j < xsec_sig.size(); ++j){
    TH2F *h = (TH2F*)file.at(0)->Get(histname)->Clone( (histname + std::to_string( j + file.size()) ).c_str() );
    hist_sig.push_back(h);
  }

  // Vector of Chi2 with different xsec_sig
  Double_t Chi2j = 0.;
  std::vector<Double_t> Chi2;

  // Histogram of sum of bkgs
  TH2F *hist_bkg = (TH2F*)hist.at(1)->Clone( (histname +  std::to_string( xsec_sig.size() + file.size()) ).c_str() );
  for(Int_t i = 2; i < hist.size(); ++i){
    hist_bkg->Add(hist.at(i));
  }

  for(Int_t j = 0; j < xsec_sig.size(); ++j){

    // Scale properly the signal
    weights.at(0) = lumi*xsec_sig.at(j) / ntot.at(0);
    hist_sig.at(j)->Sumw2();
    hist_sig.at(j)->Scale(weights.at(0));

    // Calculate Chi2
    Chi2j = 0.;

    for(Int_t kx = 1; kx <= hist_sig.at(j)->GetXaxis()->GetNbins(); ++kx){
      for(Int_t ky = 1; ky <= hist_sig.at(j)->GetYaxis()->GetNbins(); ++ky){
        Chi2j += pow( hist_sig.at(j)->GetBinContent(kx,ky) ,2) / hist_bkg->GetBinContent(kx,ky) ;
      }
    }

    Chi2.push_back(Chi2j);
  }

  // Graph of Chi2 vs xsec_sig

  TGraph *cl_gr = new TGraph();
  cl_gr->SetName(TString::Format("CLgraph_%s", histname));
  cl_gr->SetMarkerStyle(kFullCircle);
  cl_gr->SetMarkerColor(kBlack);
  cl_gr->SetLineColor(kGreen+2);
  cl_gr->SetLineWidth(2);
  cl_gr->GetXaxis()->SetTitleSize(.045);
  cl_gr->GetYaxis()->SetTitleSize(.045);
  cl_gr->GetXaxis()->SetTitle("#sigma_{sig} [pb]");
  cl_gr->GetYaxis()->SetTitle("#chi^{2}");

  for(Int_t j = 0; j < xsec_sig.size(); ++j){
    cl_gr->AddPoint(xsec_sig.at(j),Chi2.at(j));
  }

  // CL limit

  Double_t CL = TMath::ChisquareQuantile(0.95, hist_bkg->GetXaxis()->GetNbins()*hist_bkg->GetYaxis()->GetNbins());

  TLine* hor = new TLine(xsec_sig.at(0) , CL, xsec_sig.at(xsec_sig.size()-1), CL);
  hor->SetLineColor(kRed);
  hor->SetLineWidth(2);

  TString canva;
  canva.Form("cCL%s",histname);

  TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  //gPad->SetLogy();
 
  cl_gr->Draw("APL");
  hor->Draw("L same");

  // Legend
  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, Category %.1s}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L= %.0f fb^{-1}}}{ 95%% CL limits on #sigma_{sig} using " + legtitle + "}", outname, lumi/1000.);

  legend->SetHeader(head,"C");
  CustomizeLeg(legend);

  // Add the right entry to the legend
  legend->AddEntry(cl_gr,"Stat only");
  legend->AddEntry(hor,"95% CL");
  legend->Draw();
 
  TString outpdf;
  outpdf.Form(outfolder + "CL_%s.pdf",outname);
  c->SaveAs(outpdf);

  //c->Close();  // Write graph on file
  cl_gr->Write("",TObject::kSingleKey);
}
//===============================
// Function to calculate and plot Chi2 and CL limits with systematics
// File vector must have signal file first and then all bkg files
void Chi2CL2Dsys(std::vector<TFile*> file,TFile* sysfile,std::vector<TString> systag, const char *histname,std::vector<Double_t> xsec_sig,std::vector<Double_t> xsec,std::vector<Double_t> ntot,Double_t lumi,std::array<Double_t,4> leglimits,TString legtitle,const char *outname, TString outfolder){

   // Control size of inputs
   if( !( (file.size() == xsec.size()) && (file.size() == ntot.size()) )  ){
    cout << "Error in Chi2CL2Dsys: number of histograms not defined in " << histname << endl;
    cout << "file size: " << file.size() << "   xsec size: " << xsec.size() << endl;
    cout << "ntot size: " << ntot.size() << endl;
    exit(EXIT_FAILURE);
  }

  // Take histograms of 1sigma systematics
  std::vector<TH2F*> syshist;

  for(Int_t i = 0; i < systag.size(); ++i){
    TH2F *hs = (TH2F*)sysfile->Get( TString::Format( "%s_sys" + systag.at(i), histname).Data() )->Clone( TString::Format( "%s_sys" + systag.at(i), histname).Data() );
    syshist.push_back(hs);
    cout << "Got Systematic " << systag.at(i) << endl;
  }

  // Take histograms from files
  std::vector<TH2F*> hist;

  for(Int_t i = 0; i < file.size(); ++i){
    TH2F *h = (TH2F*)file.at(i)->Get(histname)->Clone( (histname + std::to_string(i)).c_str() );
    hist.push_back(h);
  }

  // Scale histograms by Integrated Luminosity, Cross-section and number of generated events
  std::vector<Double_t> weights;

  for(Int_t i = 0; i < file.size(); ++i){
    weights.push_back( lumi*xsec.at(i) / ntot.at(i) );
  }

  for(Int_t i = 1; i < hist.size(); ++i){
    hist.at(i)->Sumw2();
    hist.at(i)->Scale(weights.at(i));
  }

  // Calculate Chi2 for each value of xsec_sig

  // Vector of signal histograms with different xsec_sig
  std::vector<TH2F*> hist_sig;

  for(Int_t j = 0; j < xsec_sig.size(); ++j){
    TH2F *h = (TH2F*)file.at(0)->Get(histname)->Clone( (histname + std::to_string( j + file.size()) ).c_str() );
    hist_sig.push_back(h);
  }

  // Vector of Chi2 with different xsec_sig
  Double_t Chi2j = 0.;
  std::vector<Double_t> Chi2;
  Double_t syssum = 0.;
  Double_t sysChi2j = 0.;
  std::vector<Double_t> sysChi2;

  // Histogram of sum of bkgs
  TH2F *hist_bkg = (TH2F*)hist.at(1)->Clone( (histname +  std::to_string( xsec_sig.size() + file.size()) ).c_str() );
  for(Int_t i = 2; i < hist.size(); ++i){
    hist_bkg->Add(hist.at(i));
  }

  for(Int_t j = 0; j < xsec_sig.size(); ++j){

    // Scale properly the signal
    weights.at(0) = lumi*xsec_sig.at(j) / ntot.at(0);
    hist_sig.at(j)->Sumw2();
    hist_sig.at(j)->Scale(weights.at(0));

    // Calculate Chi2
    Chi2j = 0.;
    sysChi2j = 0.;

    for(Int_t kx = 1; kx <= hist_sig.at(j)->GetXaxis()->GetNbins(); ++kx){
      for(Int_t ky = 1; ky <= hist_sig.at(j)->GetYaxis()->GetNbins(); ++ky){

        // Stat uncertainty only
        Chi2j += pow( hist_sig.at(j)->GetBinContent(kx,ky) ,2) / hist_bkg->GetBinContent(kx,ky) ;
        
        // Stat + Sys uncertainty
        syssum = 0.;
        for(Int_t s = 0; s < systag.size(); ++s){
          syssum += pow( syshist.at(s)->GetBinContent(kx,ky) ,2) ;
        }  
        sysChi2j += pow( hist_sig.at(j)->GetBinContent(kx,ky) ,2) / (hist_bkg->GetBinContent(kx,ky) + syssum) ;
      }
    }

    Chi2.push_back(Chi2j);
    sysChi2.push_back(sysChi2j);
  }

  // Graphs of Chi2 vs xsec_sig
  // Stat only
  TGraph *cl_gr = new TGraph();
  cl_gr->SetName(TString::Format("CLgraph_%s", histname));
  cl_gr->SetMarkerStyle(kFullCircle);
  cl_gr->SetMarkerColor(kBlack);
  cl_gr->SetLineColor(kGreen+2);
  cl_gr->SetLineWidth(2);
  cl_gr->GetXaxis()->SetTitleSize(.045);
  cl_gr->GetYaxis()->SetTitleSize(.045);
  cl_gr->GetXaxis()->SetTitle("#sigma_{sig} [pb]");
  cl_gr->GetYaxis()->SetTitle("#chi^{2}");

  // Stat + Sys
  TGraph *syscl_gr = new TGraph();
  syscl_gr->SetName(TString::Format("sysCLgraph_%s", histname));
  syscl_gr->SetMarkerStyle(kFullCircle);
  syscl_gr->SetMarkerColor(kBlack);
  syscl_gr->SetLineColor(kViolet-1);
  syscl_gr->SetLineWidth(2);
  syscl_gr->GetXaxis()->SetTitleSize(.045);
  syscl_gr->GetYaxis()->SetTitleSize(.045);
  syscl_gr->GetXaxis()->SetTitle("#sigma_{sig} [pb]");
  syscl_gr->GetYaxis()->SetTitle("#chi^{2}");

  // Fill graphs
  for(Int_t j = 0; j < xsec_sig.size(); ++j){
    cl_gr->AddPoint(xsec_sig.at(j),Chi2.at(j));
    syscl_gr->AddPoint(xsec_sig.at(j),sysChi2.at(j));
  }

  // CL limit

  Double_t CL = TMath::ChisquareQuantile(0.95, hist_bkg->GetXaxis()->GetNbins()*hist_bkg->GetYaxis()->GetNbins());

  TLine* hor = new TLine(xsec_sig.at(0) , CL, xsec_sig.at(xsec_sig.size()-1), CL);
  hor->SetLineColor(kRed);
  hor->SetLineWidth(2);

  TString canva;
  canva.Form("cCLsys%s",histname);

  TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleXSize(0.045);
  gStyle->SetTitleYSize(0.045);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();

  TMultiGraph* mg = new TMultiGraph;
  mg->SetTitle(";#sigma_{sig} [pb];#chi^{2}");
  mg->Add(cl_gr);
  mg->Add(syscl_gr);
  mg->Draw("APL"); 

  //cl_gr->Draw("APL");
  //syscl_gr->Draw("APL");

  hor->Draw("L same");

  // Legend
  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, Category %.1s}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L= %.0f fb^{-1}}}{ 95%% CL limits on #sigma_{sig} using " + legtitle + "}", outname, lumi/1000.);

  legend->SetHeader(head,"C");
  CustomizeLeg(legend);

  // Add the right entry to the legend
  legend->AddEntry(cl_gr,"Stat only");
  legend->AddEntry(syscl_gr,"Stat + Sys");
  legend->AddEntry(hor,"95% CL");
  legend->Draw();
 
  TString outpdf;
  outpdf.Form(outfolder + "sysCL_%s.pdf",outname);
  c->SaveAs(outpdf);

  //c->Close();  
  
  // Write graph on file
  cl_gr->Write("",TObject::kSingleKey);
  syscl_gr->Write("",TObject::kSingleKey);
}
//===============================
// Function to calculate and plot Chi2 and CL limits with single and combined systematics
// File vector must have signal file first and then all bkg files
void Chi2CL2Dsinglesys(std::vector<TFile*> file,TFile* sysfile,std::vector<TString> systag, const char *histname,std::vector<Double_t> xsec_sig,std::vector<Double_t> xsec,std::vector<Double_t> ntot,Double_t lumi,std::array<Double_t,4> leglimits,TString legtitle,const char *outname, TString outfolder, std::vector<Color_t> singlesyscolors,std::vector<TString> syslegentry){

  // Control size of inputs
  if( !( (file.size() == xsec.size()) && (file.size() == ntot.size()) )  ){
    cout << "Error in Chi2CL2Dsinglesys: number of histograms not defined in " << histname << endl;
    cout << "file size: " << file.size() << "   xsec size: " << xsec.size() << endl;
    cout << "ntot size: " << ntot.size() << endl;
    exit(EXIT_FAILURE);
  }
  if( systag.size() != singlesyscolors.size() ){
    cout << "Error in Chi2CL2Dsinglesys: number of systematics not defined in " << histname << endl;
    cout << "systag size: " << systag.size() << "   singlesyscolors size: " << singlesyscolors.size() << endl;
    exit(EXIT_FAILURE);
  }

  // Take histograms of 1sigma systematics
  std::vector<TH2F*> syshist;

  for(Int_t i = 0; i < systag.size(); ++i){
    TH2F *hs = (TH2F*)sysfile->Get( TString::Format( "%s_sys" + systag.at(i), histname).Data() )->Clone( TString::Format( "%s_sys" + systag.at(i), histname).Data() );
    syshist.push_back(hs);
    cout << "Got Systematic " << systag.at(i) << endl;
  }

  // Take histograms from files
  std::vector<TH2F*> hist;

  for(Int_t i = 0; i < file.size(); ++i){
    TH2F *h = (TH2F*)file.at(i)->Get(histname)->Clone( (histname + std::to_string(i)).c_str() );
    hist.push_back(h);
  }

  // Scale histograms by Integrated Luminosity, Cross-section and number of generated events
  std::vector<Double_t> weights;

  for(Int_t i = 0; i < file.size(); ++i){
    weights.push_back( lumi*xsec.at(i) / ntot.at(i) );
  }

  for(Int_t i = 1; i < hist.size(); ++i){
    hist.at(i)->Sumw2();
    hist.at(i)->Scale(weights.at(i));
  }

  // Calculate Chi2 for each value of xsec_sig

  // Vector of signal histograms with different xsec_sig
  std::vector<TH2F*> hist_sig;

  for(Int_t j = 0; j < xsec_sig.size(); ++j){
    TH2F *h = (TH2F*)file.at(0)->Get(histname)->Clone( (histname + std::to_string( j + file.size()) ).c_str() );
    hist_sig.push_back(h);
  }

  // Vector of Chi2 with different xsec_sig
  Double_t Chi2j = 0.;
  std::vector<Double_t> Chi2;
  Double_t syssum = 0.;
  Double_t sysChi2j = 0.;
  std::vector<Double_t> sysChi2;
  std::vector<Double_t> singlesysChi2j;
  std::vector<std::vector<Double_t>> singlesysChi2;
  for(Int_t s = 0; s < systag.size(); ++s){
    singlesysChi2j.push_back(0.);
    singlesysChi2.push_back({});
  }

  // Histogram of sum of bkgs
  TH2F *hist_bkg = (TH2F*)hist.at(1)->Clone( (histname +  std::to_string( xsec_sig.size() + file.size()) ).c_str() );
  for(Int_t i = 2; i < hist.size(); ++i){
    hist_bkg->Add(hist.at(i));
  }

  for(Int_t j = 0; j < xsec_sig.size(); ++j){

    // Scale properly the signal
    weights.at(0) = lumi*xsec_sig.at(j) / ntot.at(0);
    hist_sig.at(j)->Sumw2();
    hist_sig.at(j)->Scale(weights.at(0));

    // Calculate Chi2
    Chi2j = 0.;
    sysChi2j = 0.;
    for(Int_t c = 0; c < systag.size(); ++c){
      singlesysChi2j.at(c) = 0.;
    }

    for(Int_t kx = 1; kx <= hist_sig.at(j)->GetXaxis()->GetNbins(); ++kx){
      for(Int_t ky = 1; ky <= hist_sig.at(j)->GetYaxis()->GetNbins(); ++ky){

        // Stat uncertainty only
        Chi2j += pow( hist_sig.at(j)->GetBinContent(kx,ky) ,2) / hist_bkg->GetBinContent(kx,ky) ;

        // Stat + Single Sys uncertainty
        for(Int_t s = 0; s < systag.size(); ++s){
          singlesysChi2j.at(s) += pow( hist_sig.at(j)->GetBinContent(kx,ky) ,2) / (hist_bkg->GetBinContent(kx,ky) + pow( syshist.at(s)->GetBinContent(kx,ky) ,2)) ;
        }  
        
        // Stat + Sys uncertainty
        syssum = 0.;
        for(Int_t s = 0; s < systag.size(); ++s){
          syssum += pow( syshist.at(s)->GetBinContent(kx,ky) ,2) ;
        }  
        sysChi2j += pow( hist_sig.at(j)->GetBinContent(kx,ky) ,2) / (hist_bkg->GetBinContent(kx,ky) + syssum) ;
      }
    }

    Chi2.push_back(Chi2j);
    sysChi2.push_back(sysChi2j);
    for(Int_t s = 0; s < systag.size(); ++s){
      singlesysChi2.at(s).push_back(singlesysChi2j.at(s));
    }
  }

  // Graphs of Chi2 vs xsec_sig
  // Stat only
  TGraph *cl_gr = new TGraph();
  cl_gr->SetName(TString::Format("CLgraph_%s", histname));
  CustomizeGr(cl_gr, kGreen+2);

  // Stat + Single Sys
  std::vector<TGraph *> singlesyscl_gr;

  for(Int_t s = 0; s < systag.size(); ++s){
    TGraph *gr = new TGraph();
    gr->SetName(TString::Format("singlesysCLgraph_%s" + systag.at(s), histname));
    CustomizeGr(gr, singlesyscolors.at(s));
    singlesyscl_gr.push_back(gr);
  }

  // Stat + Sys
  TGraph *syscl_gr = new TGraph();
  syscl_gr->SetName(TString::Format("sysCLgraph_%s", histname));
  CustomizeGr(syscl_gr, kViolet-1);

  // Fill graphs
  for(Int_t j = 0; j < xsec_sig.size(); ++j){
    cl_gr->AddPoint(xsec_sig.at(j),Chi2.at(j));
    syscl_gr->AddPoint(xsec_sig.at(j),sysChi2.at(j));
    for(Int_t s = 0; s < systag.size(); ++s){
      singlesyscl_gr.at(s)->AddPoint(xsec_sig.at(j),singlesysChi2.at(s).at(j));
    }
  }

  // CL limit

  Double_t CL = TMath::ChisquareQuantile(0.95, hist_bkg->GetXaxis()->GetNbins()*hist_bkg->GetYaxis()->GetNbins());

  TLine* hor = new TLine(xsec_sig.at(0) , CL, xsec_sig.at(xsec_sig.size()-1), CL);
  hor->SetLineColor(kRed);
  hor->SetLineWidth(2);

  TString canva;
  canva.Form("cCLsinglesys%s",histname);

  TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleXSize(0.045);
  gStyle->SetTitleYSize(0.045);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();

  TMultiGraph* mg = new TMultiGraph;
  mg->SetTitle(";#sigma_{sig} [pb];#chi^{2}");
  mg->Add(cl_gr);
  mg->Add(syscl_gr);
  for(Int_t s = 0; s < systag.size(); ++s){
    mg->Add(singlesyscl_gr.at(s));
  }
  mg->Draw("APL"); 

  hor->Draw("L same");

  // Legend
  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, Category %.1s}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L= %.0f fb^{-1}}}{ 95%% CL limits on #sigma_{sig} using " + legtitle + "}", outname, lumi/1000.);

  legend->SetHeader(head,"C");
  CustomizeLeg(legend);
  legend->SetNColumns(2);

  // Add the right entry to the legend
  legend->AddEntry(cl_gr,"Stat only");
  legend->AddEntry(syscl_gr,"Stat + Sys");
  for(Int_t s = 0; s < systag.size(); ++s){
    legend->AddEntry(singlesyscl_gr.at(s),"Stat + " + syslegentry.at(s));
  }
  legend->AddEntry(hor,"95% CL");
  legend->Draw();
 
  TString outpdf;
  outpdf.Form(outfolder + "singlesysCL_%s.pdf",outname);
  c->SaveAs(outpdf);

  //c->Close();  
  
  // Write graph on file
  cl_gr->Write("",TObject::kSingleKey);
  syscl_gr->Write("",TObject::kSingleKey);
  for(Int_t s = 0; s < systag.size(); ++s){
    singlesyscl_gr.at(s)->Write("",TObject::kSingleKey);
  }
}
//===============================
//---------- END FUNCTIONS -------------------


void Chi2_MergeHistos()
{

  gSystem->Load("libDelphes");

  // Import histograms 
  TString dir = "/gfsvol01/atlas/giuliac/plots_and_outputs/20250718_Task27/";
  TString sampledir = dir;
  TString outdir = "/gfsvol01/atlas/giuliac/plots_and_outputs/20250726_Task28/";

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

  // File with 1D systematics histograms
  TString sysname = sampledir + "sys.root";
  TFile *f_sys = new TFile(sysname);

  // Systematics tags
  std::vector<TString> systag = {"_JES", "_etrk", "_eb", "_xsec"};
  std::vector<TString> syslegentry = {"Jet Energy Scale", "Charged Hadron Tracking Efficiency", "b-tagging Efficiency", "Z+jets Background Estimation"};

  // Integrated Luminosity [pb^-1] (arbitrary, no real data here)
  Double_t lumi = 360000;

  // Cross-section [pb] 
  //Double_t xsec_bkgZq = 1.250e4;
  //Double_t xsec_bkgZg = 7.759e3;
  Double_t xsec_bkgZq = 1.671e3;
  Double_t xsec_bkgZg = 4.374e2;
  Double_t xsec_bkgZZ = 7.688e-1;
  Double_t xsec_bkgHZ = 6.455e-2;
  Double_t xsec_sig = xsec_bkgZq + xsec_bkgZg +xsec_bkgZZ + xsec_bkgHZ;

  // Number of MC events
  Double_t ntot_sig = 500000;
  Double_t ntot_bkgZq = 5000000;
  Double_t ntot_bkgZg = 5000000;
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
  sampledir = outdir;
  sampledir += "area_" ;

  // write all scaled histograms on file
  TString outname = sampledir + "MergedHistos.root";

  // Event categories
  std::vector<TString> cat = {"A_","B_"};

  TFile *f = new TFile(outname,"RECREATE");

  std::array<Double_t,4> leglimits_1 = {0.461153,0.548611,0.77193,0.828125};
  std::array<Double_t,4> leglimits2D_1 = {0.190476,0.668403,0.839599,0.878472};
  std::array<Double_t,4> leglimits2D_sys = {0.195489,0.755208,0.824561,0.885417};

  Bool_t area = kTRUE;

  for(Int_t i = 0; i < cat.size(); ++i){

    /*
    // 1D histograms

    DrawHistos(area,infiles,cat.at(i) + "Invariant_mass_of_leading_and_subleading_small_jets",weights,colors,"m_{jj} [GeV/c^{2}]","#events",leglimits_1,"m_{jj} of first 2 Small Jets",legentry,cat.at(i) + "InvMass_small",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "Invariant_mass_of_leading_and_subleading_large_jets",weights,colors,"m_{jj} [GeV/c^{2}]","#events",leglimits_1,"m_{jj} of first 2 Large Jets",legentry,cat.at(i) + "InvMass_large",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "Invariant_mass_of_all_btagged_small_jets",weights,colors,"Invariant Mass [GeV/c^{2}]","#events",leglimits_1,"Invariant Mass of all b-tagged Small Jets",legentry,cat.at(i) + "InvMass_allb",sampledir);

    DrawHistos(area,infiles,cat.at(i) + "DeltaR_between_leading_and_subleading_small_jets",weights,colors,"#Delta R","# events",leglimits_1,"#Delta R between first 2 Small Jets",legentry,cat.at(i) + "DeltaRS",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "DeltaPhi_between_leading_and_subleading_small_jets",weights,colors,"|#Delta #varphi|","# events",leglimits_1,"|#Delta #varphi| between first 2 Small Jets",legentry,cat.at(i) + "DeltaPhiS",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "DeltaR_between_leading_and_subleading_large_jets",weights,colors,"#Delta R","# events",leglimits_1,"#Delta R between first 2 Large Jets",legentry,cat.at(i) + "DeltaRL",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "DeltaPhi_between_leading_and_subleading_large_jets",weights,colors,"|#Delta #varphi|","# events",leglimits_1,"|#Delta #varphi| between first 2 Large Jets",legentry,cat.at(i) + "DeltaPhiL",sampledir);

    DrawHistos(area,infiles,cat.at(i) + "leading_small_jet_pt",weights,colors,"p_{T} [GeV/c]","# events",leglimits_1,"pT of Leading Small Jet",legentry,cat.at(i) + "JetPT1S",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "subleading_small_jet_pt",weights,colors,"p_{T} [GeV/c]","# events",leglimits_1,"pT of Subleading Small Jet",legentry,cat.at(i) + "JetPT2S",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "leading_large_jet_pt",weights,colors,"p_{T} [GeV/c]","# events",leglimits_1,"pT of Leading Large Jet",legentry,cat.at(i) + "JetPT1L",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "subleading_large_jet_pt",weights,colors,"p_{T} [GeV/c]","# events",leglimits_1,"pT of Subleading Large Jet",legentry,cat.at(i) + "JetPT2L",sampledir);

    DrawHistos(area,infiles,cat.at(i) + "number_of_constituents_of_leading_large_jet",weights,colors,"Number of Constituents","# events",leglimits_1,"Number of Constituents of Leading Large Jet",legentry,cat.at(i) + "NConst1",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "number_of_constituents_of_subleading_large_jet",weights,colors,"Number of Constituents","# events",leglimits_1,"Number of Constituents of Subleading Large Jet",legentry,cat.at(i) + "NConst2",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "number_of_small_jets",weights,colors,"Number of Jets","# events",leglimits_1,"Number of Small Jets",legentry,cat.at(i) + "Nsmall",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "number_of_btagged_small_jets",weights,colors,"Number of Jets","# events",leglimits_1,"Number of b-tagged Small Jets",legentry,cat.at(i) + "Nbsmall",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "number_of_large_jets",weights,colors,"Number of Jets","# events",leglimits_1,"Number of Large Jets",legentry,cat.at(i) + "Nlarge",sampledir);

    DrawHistos(area,infiles,cat.at(i) + "Prompt_Track_Fraction_of_leading_large_jet",weights,colors,"PTF","# events",leglimits_1,"PTF of Leading Large Jet",legentry,cat.at(i) + "PTF1L",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "Prompt_Track_Fraction_of_subleading_large_jet",weights,colors,"PTF","# events",leglimits_1,"PTF of Subleading Large Jet",legentry,cat.at(i) + "PTF2L",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "Prompt_Track_Fraction_of_leading_small_jet",weights,colors,"PTF","# events",leglimits_1,"PTF of Leading Small Jet",legentry,cat.at(i) + "PTF1S",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "Prompt_Track_Fraction_of_subleading_small_jet",weights,colors,"PTF","# events",leglimits_1,"PTF of Subleading Small Jet",legentry,cat.at(i) + "PTF2S",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "Prompt_Track_Fraction_of_leading_btagged_small_jet",weights,colors,"PTF","# events",leglimits_1,"PTF of Leading b-tagged Small Jet",legentry,cat.at(i) + "PTF1bS",sampledir);

    DrawHistos(area,infiles,cat.at(i) + "2point_Energy_Correlation_div_pT_of_leading_large_jet",weights,colors,"ECF2/pT [GeV/c]","# events",leglimits_1,"ECF2/pT of Leading Large Jet",legentry,cat.at(i) + "ECF1L",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "2point_Energy_Correlation_div_pT_of_subleading_large_jet",weights,colors,"ECF2/pT [GeV/c]","# events",leglimits_1,"ECF2/pT of Subleading Large Jet",legentry,cat.at(i) + "ECF2L",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "2point_Energy_Correlation_div_pT_of_leading_small_jet",weights,colors,"ECF2/pT [GeV/c]","# events",leglimits_1,"ECF2/pT of Leading Small Jet",legentry,cat.at(i) + "ECF1S",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "2point_Energy_Correlation_div_pT_of_subleading_small_jet",weights,colors,"ECF2/pT [GeV/c]","# events",leglimits_1,"ECF2/pT of Subleading Small Jet",legentry,cat.at(i) + "ECF2S",sampledir);

    DrawHistos(area,infiles,cat.at(i) + "pt_of_first_two_leptons",weights,colors,"p_{Tll} [GeV/c]","# events",leglimits_1,"pT_{ll} of Z tagging Leptons",legentry,cat.at(i) + "Lep2PT",sampledir);
    DrawHistos(area,infiles,cat.at(i) + "DeltaR_between_first_two_leptons",weights,colors,"#Delta R","# events",leglimits_1,"#Delta R between Z tagging Leptons",legentry,cat.at(i) + "LepDeltaR",sampledir);

    // 2D histograms

    DrawHistos2D(area,infiles,cat.at(i) + "IMsmall_vs_DeltaRsmall",weights,"m_{jj} of first 2 Small Jets [GeV/c^2]","#Delta R between first 2 Small Jets","# events",leglimits2D_1,"m_{jj} vs #Delta R of first 2 Small Jets",legentry,cat.at(i) + "IMs_DRs",sampledir);
    DrawHistos2D(area,infiles,cat.at(i) + "IMlarge_vs_DeltaRlarge",weights,"m_{jj} of first 2 Large Jets [GeV/c^2]","#Delta R between first 2 Large Jets","# events",leglimits2D_1,"m_{jj} vs #Delta R of first 2 Large Jets",legentry,cat.at(i) + "IMl_DRl",sampledir);

    DrawHistos2D(area,infiles,cat.at(i) + "PTF1small_vs_PTF2small",weights,"PTF of Leading Small Jet","PTF of Subleading Small Jet","# events",leglimits2D_1,"Prompt Track Fraction of Leading vs Subleading Small Jet",legentry,cat.at(i) + "PTF1S_PTF2S",sampledir);
    DrawHistos2D(area,infiles,cat.at(i) + "PTF1small_vs_PTF1large",weights,"PTF of Leading Small Jet","PTF of Leading Large Jet","# events",leglimits2D_1,"Prompt Track Fraction of Leading Small vs Large Jet",legentry,cat.at(i) + "PTF1S_PTF1L",sampledir);
    DrawHistos2D(area,infiles,cat.at(i) + "PTF1bsmall_vs_PTF1large",weights,"PTF of Leading b-tagged Small Jet","PTF of Leading Large Jet","# events",leglimits2D_1,"Prompt Track Fraction of Leading b-tagged Small vs Leading Large Jet",legentry,cat.at(i) + "PTF1bS_PTF1L",sampledir);
    
    // Histos for Chi2 calculation
    DrawHistos(area,infiles,"Chi2" + cat.at(i) + "Invariant_mass_of_leading_and_subleading_small_jets",weights,colors,"m_{jj} [GeV/c^{2}]","#events",leglimits_1,"m_{jj} of first 2 Small Jets",legentry,cat.at(i) + "InvMass_small" + "_Chi2",sampledir);
    DrawHistos2D(area,infiles,"Chi2HD" + cat.at(i) + "PTF1small_vs_PTF2small",weights,"PTF of Leading Small Jet","PTF of Subleading Small Jet","# events",leglimits2D_1,"Prompt Track Fraction of Leading vs Subleading Small Jet",legentry,cat.at(i) + "PTF1S_PTF2S" + "_Chi2HD",sampledir);
    DrawHistos2D(area,infiles,"Chi2LD" + cat.at(i) + "PTF1small_vs_PTF2small",weights,"PTF of Leading Small Jet","PTF of Subleading Small Jet","# events",leglimits2D_1,"Prompt Track Fraction of Leading vs Subleading Small Jet",legentry,cat.at(i) + "PTF1S_PTF2S" + "_Chi2LD",sampledir);
    */

    // Histos of 1sigma systematics
    for(Int_t j = 0; j < systag.size(); ++j){
      DrawSys(f_sys, "Chi2LD" + cat.at(i) + "PTF1small_vs_PTF2small" + "_sys" + systag.at(j),"PTF of Leading Small Jet","PTF of Subleading Small Jet","# events",leglimits2D_sys,"1#sigma deviation on total background due to systematic",syslegentry.at(j),cat.at(i) + "sys" + systag.at(j),outdir);
    }

  }


  // Calculation of Chi2 and CL limits

  // Signal xsections [pb]
  std::vector<Double_t> xsec_s1;
  for(Double_t i = 0.08; i < 0.122; i += 0.002){
    xsec_s1.push_back(i);
  }

  std::vector<Double_t> xsec_s2;
  for(Double_t i = 0.04; i < 0.62; i += 0.02){
    xsec_s2.push_back(i);
  }

  std::array<Double_t,4> leglimits_2 = {0.294486,0.664931,0.605263,0.835069};
  std::array<Double_t,4> leglimits_3 = {0.293233,0.604167,0.60401,0.835069};
  std::array<Double_t,4> leglimits_4 = {0.484962,0.142361,0.803258,0.369792};
  std::array<Double_t,4> leglimits_5 = {0.484962,0.128472,0.803258,0.373264};

  std::vector<Color_t> singlesyscolors = {kCyan-6,kOrange-3,kAzure+7,kMagenta-6};
  std::vector<TString> shortsyslegentry = {"JES","#varepsilon_{trk}","#varepsilon_{b}","#sigma_{Z+jets}"};

  for(Int_t i = 0; i < cat.size(); ++i){
    Chi2CL2D(infiles,"Chi2LD" + cat.at(i) + "PTF1small_vs_PTF2small",xsec_s1,xsec,ntot,lumi,leglimits_2,"PTF of L vs Subl Small Jet",cat.at(i) + "PTF1S_PTF2S" + "_LD", outdir);
    Chi2CL2Dsys(infiles,f_sys,systag,"Chi2LD" + cat.at(i) + "PTF1small_vs_PTF2small",xsec_s2,xsec,ntot,lumi,leglimits_4,"PTF of L vs Subl Small Jet",cat.at(i) + "PTF1S_PTF2S" + "_LD",outdir);
    Chi2CL2Dsinglesys(infiles,f_sys,systag,"Chi2LD" + cat.at(i) + "PTF1small_vs_PTF2small",xsec_s2,xsec,ntot,lumi,leglimits_5,"PTF of L vs Subl Small Jet",cat.at(i) + "PTF1S_PTF2S" + "_LD",outdir,singlesyscolors,shortsyslegentry);
  }


}//end MergeHistos
