/*
Author: Giulia Cossutti

Macro to combine calculation of Chi2 and 95% CL limits for cathegories A and B
with different signal xsections.
Both Stat only and Stat+Sys combinations are calculated.

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
root -l /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250718_Task27/CombineCL.C 

*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "/gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/classes/DelphesClasses.h"
#include "/gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#endif

#include <numeric>

//------------------------------------------------------------------------------
//---------- FUNCTIONS -------------------
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
// Function to combine CL limits (stat only)
void CombineCLstat(TFile *infile,std::vector<TString> namecycle,std::vector<TString> cat,std::vector<TString> histtype,std::vector<TString> histname,std::vector<TString> histprefix,std::vector<TString> graphprefix,std::array<Double_t,4> leglimits,Double_t lumi,TString legtitle,TString outfolder){

  // Size control on vectors
  if( !((cat.size() == histtype.size()) && (histname.size() == histprefix.size()) && (cat.size()== graphprefix.size()) && (cat.size()==histname.size()) && (namecycle.size()== graphprefix.size())) ){
    cout << "Error in CombineCLstat: number of histograms not defined" << endl;
    cout << "cat.size(): " << cat.size() << "  histtype.size(): " << histtype.size() << "  namecycle.size(): " << namecycle.size() << endl;
    cout << "histname.size(): " << histname.size() << "  histprefix.size(): " << histprefix.size() << "  graphprefix.size(): " << graphprefix.size() << endl;
    exit(EXIT_FAILURE);
  }

  std::vector<Int_t> ndf;
  std::vector<TGraph*> cl_gr;

  for(Int_t i = 0; i < cat.size(); ++i){

    // Get degrees of freedom for the histograms
    if(std::strcmp(histtype.at(i),"TH1F") == 0){
      TH1F *h = (TH1F*)infile->Get(histprefix.at(i) + cat.at(i) + histname.at(i) + "0")->Clone( cat.at(i) + "h" );
      ndf.push_back(h->GetXaxis()->GetNbins());
      h->Write("",TObject::kSingleKey);
    }else if(std::strcmp(histtype.at(i),"TH2F") == 0){
      TH2F *h = (TH2F*)infile->Get(histprefix.at(i) + cat.at(i) + histname.at(i) + "0")->Clone( cat.at(i) + "h" );
      ndf.push_back(h->GetXaxis()->GetNbins()*h->GetYaxis()->GetNbins());
      h->Write("",TObject::kSingleKey);
    }

    // Get xsec_sig vs Chi2 graphs
    TGraph *gr = (TGraph*)infile->Get(graphprefix.at(i) + cat.at(i) + histname.at(i) + namecycle.at(i))->Clone( cat.at(i) + "grstatonly" );
    cl_gr.push_back(gr);
    gr->Write("",TObject::kSingleKey);
  }

  // Build combination graph
  TGraph *comb_gr = (TGraph*)cl_gr.at(0)->Clone( "comb_gr" );
  Double_t newChi2 = 0.;

  for(Int_t i = 1; i < cat.size(); ++i){
    // Control that graphs have same number of points
    if(comb_gr->GetN() != cl_gr.at(i)->GetN()){
      cout << "Error in CombineCLstat: different number of points in graphs" << endl;
      exit(EXIT_FAILURE);
    }

    // Add graphs point per point
    for(Int_t k = 0; k < comb_gr->GetN(); ++k){
      // Control to use same set of xsec_sig
      if( abs(comb_gr->GetPointX(k) - cl_gr.at(i)->GetPointX(k)) > 1e-5){
        cout << "Error in CombineCLstat: different PointX in graphs" << endl;
        exit(EXIT_FAILURE);
      }

      newChi2 = comb_gr->GetPointY(k) + cl_gr.at(i)->GetPointY(k);
      comb_gr->SetPointY(k,newChi2);
    }
  }

  // CL limit
  Double_t CL = TMath::ChisquareQuantile(0.95, std::accumulate(ndf.begin(), ndf.end(), 0));

  TLine* hor = new TLine(comb_gr->GetPointX(0) , CL, comb_gr->GetPointX(comb_gr->GetN()-1), CL);
  hor->SetLineColor(kRed);
  hor->SetLineWidth(2);

  // Draw combined graph
  TString canva = "cCL_combined";

  TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  //gPad->SetLogy();
 
  comb_gr->Draw("APL");
  hor->Draw("L same");

  comb_gr->Write("",TObject::kSingleKey);


  // Legend

  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, Combined Categories}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L= %.0f fb^{-1}}}{ 95%% CL limits on #sigma_{sig} using " + legtitle + "}", lumi/1000.);

  legend->SetHeader(head,"C");
  CustomizeLeg(legend);

  // Add the right entry to the legend
  legend->AddEntry(comb_gr,"Stat only");
  legend->AddEntry(hor,"95% CL");
  legend->Draw();
 
  TString outpdf = outfolder + "CombinedCLstat.pdf";
  c->SaveAs(outpdf);

  //c->Close();

}
//===================================
// Function to combine CL limits (stat+sys)
void CombineCLstatsys(TFile *infile,std::vector<TString> namecycle,std::vector<TString> sysnamecycle,std::vector<TString> cat,std::vector<TString> histtype,std::vector<TString> histname,std::vector<TString> histprefix,std::vector<TString> graphprefix,std::vector<TString> sysgraphprefix,std::array<Double_t,4> leglimits,Double_t lumi,TString legtitle,TString outfolder){

  // Size control on vectors
  if( !((cat.size() == histtype.size()) && (histname.size() == histprefix.size()) && (cat.size()== graphprefix.size()) && (cat.size()==histname.size()) && (namecycle.size()== graphprefix.size())) ){
    cout << "Error in CombineCLstatsys: number of histograms not defined" << endl;
    cout << "cat.size(): " << cat.size() << "  histtype.size(): " << histtype.size() << "  namecycle.size(): " << namecycle.size() << endl;
    cout << "histname.size(): " << histname.size() << "  histprefix.size(): " << histprefix.size() << "  graphprefix.size(): " << graphprefix.size() << endl;
    exit(EXIT_FAILURE);
  }

  std::vector<Int_t> ndf;
  std::vector<TGraph*> cl_gr;
  std::vector<TGraph*> syscl_gr;

  for(Int_t i = 0; i < cat.size(); ++i){

    // Get degrees of freedom for the histograms
    if(std::strcmp(histtype.at(i),"TH1F") == 0){
      TH1F *h = (TH1F*)infile->Get(histprefix.at(i) + cat.at(i) + histname.at(i) + "0")->Clone( cat.at(i) + "hsys" );
      ndf.push_back(h->GetXaxis()->GetNbins());
      h->Write("",TObject::kSingleKey);
    }else if(std::strcmp(histtype.at(i),"TH2F") == 0){
      TH2F *h = (TH2F*)infile->Get(histprefix.at(i) + cat.at(i) + histname.at(i) + "0")->Clone( cat.at(i) + "hsys" );
      ndf.push_back(h->GetXaxis()->GetNbins()*h->GetYaxis()->GetNbins());
      h->Write("",TObject::kSingleKey);
    }

    // Get xsec_sig vs Chi2 graphs
    TGraph *gr = (TGraph*)infile->Get(graphprefix.at(i) + cat.at(i) + histname.at(i) + namecycle.at(i))->Clone( cat.at(i) + "grstat" );
    cl_gr.push_back(gr);
    gr->Write("",TObject::kSingleKey);

    TGraph *sysgr = (TGraph*)infile->Get(sysgraphprefix.at(i) + cat.at(i) + histname.at(i) + sysnamecycle.at(i))->Clone( cat.at(i) + "grstatsys" );
    syscl_gr.push_back(sysgr);
    sysgr->Write("",TObject::kSingleKey);
  }

  // Build combination graph
  TGraph *comb_gr = (TGraph*)cl_gr.at(0)->Clone( "comb_gr" );
  Double_t newChi2 = 0.;
  TGraph *syscomb_gr = (TGraph*)syscl_gr.at(0)->Clone( "syscomb_gr" );
  Double_t sysnewChi2 = 0.;

  for(Int_t i = 1; i < cat.size(); ++i){
    // Control that graphs have same number of points
    if( (comb_gr->GetN() != cl_gr.at(i)->GetN()) || (syscomb_gr->GetN() != syscl_gr.at(i)->GetN()) || (comb_gr->GetN() != syscomb_gr->GetN())){
      cout << "Error in CombineCLstatsys: different number of points in graphs" << endl;
      exit(EXIT_FAILURE);
    }

    // Add graphs point per point
    for(Int_t k = 0; k < comb_gr->GetN(); ++k){
      // Control to use same set of xsec_sig
      if( (abs(comb_gr->GetPointX(k) - cl_gr.at(i)->GetPointX(k)) > 1e-5)  || (abs(syscomb_gr->GetPointX(k) - syscl_gr.at(i)->GetPointX(k)) > 1e-5) ){
        cout << "Error in CombineCLstatsys: different PointX in graphs" << endl;
        exit(EXIT_FAILURE);
      }

      newChi2 = comb_gr->GetPointY(k) + cl_gr.at(i)->GetPointY(k);
      sysnewChi2 = syscomb_gr->GetPointY(k) + syscl_gr.at(i)->GetPointY(k);
      comb_gr->SetPointY(k,newChi2);
      syscomb_gr->SetPointY(k,sysnewChi2);
    }
  }

  // CL limit
  Double_t CL = TMath::ChisquareQuantile(0.95, std::accumulate(ndf.begin(), ndf.end(), 0));

  TLine* hor = new TLine(comb_gr->GetPointX(0) , CL, comb_gr->GetPointX(comb_gr->GetN()-1), CL);
  hor->SetLineColor(kRed);
  hor->SetLineWidth(2);

  // Draw combined graph
  TString canva = "cCLsys_combined";

  TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleXSize(0.045);
  gStyle->SetTitleYSize(0.045);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();
 
  TMultiGraph* mg = new TMultiGraph;
  mg->SetTitle(";#sigma_{sig} [pb];#chi^{2}");
  mg->Add(comb_gr);
  mg->Add(syscomb_gr);
  mg->Draw("APL"); 

  //comb_gr->Draw("APL");
  //syscomb_gr->Draw("APL");

  hor->Draw("L same");

  comb_gr->Write("",TObject::kSingleKey);
  syscomb_gr->Write("",TObject::kSingleKey);


  // Legend

  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = TString::Format("#splitline{#splitline{Model D, c#tau_{#pi_{D}} = 50 mm, Combined Categories}{PYTHIA8+DELPHES, #sqrt{s}=13.6 TeV, L= %.0f fb^{-1}}}{ 95%% CL limits on #sigma_{sig} using " + legtitle + "}", lumi/1000.);

  legend->SetHeader(head,"C");
  CustomizeLeg(legend);

  // Add the right entry to the legend
  legend->AddEntry(comb_gr,"Stat only");
  legend->AddEntry(syscomb_gr,"Stat + Sys");
  legend->AddEntry(hor,"95% CL");
  legend->Draw();
 
  TString outpdf = outfolder + "CombinedCLstatsys.pdf";
  c->SaveAs(outpdf);

  //c->Close();

}
//===================================
//---------- END FUNCTIONS -------------------


void CombineCL()
{

  gSystem->Load("libDelphes");

  // Import histograms and graphs
  TString infilename = "/gfsvol01/atlas/giuliac/plots_and_outputs/20250718_Task27/area_MergedHistos.root";
  TFile *infile = new TFile(infilename);

  // Names of histograms from different cathegories
  std::vector<TString> cat = {"A_","B_"};
  std::vector<TString> histtype = {"TH2F","TH2F"};
  std::vector<TString> histname = {"PTF1small_vs_PTF2small","PTF1small_vs_PTF2small"};
  std::vector<TString> histprefix = {"Chi2LD","Chi2LD"};
  std::vector<TString> graphprefix = {"CLgraph_Chi2LD","CLgraph_Chi2LD"};
  std::vector<TString> namecycle = {";1",";1"};

  // Legend features
  std::array<Double_t,4> leglimits = {0.293233,0.633681,0.60401,0.868056};

  Double_t lumi = 360000;
  TString legtitle = "PTF of L vs Subl Small Jet (A,B)";

  // Outputs
  TString outfolder = "/gfsvol01/atlas/giuliac/plots_and_outputs/20250718_Task27/";
  TFile *f = new TFile(outfolder + "CL_graphs.root","RECREATE");

  // Combine Stat only limits
  CombineCLstat(infile,namecycle,cat,histtype,histname,histprefix,graphprefix,leglimits,lumi,legtitle,outfolder);

  // Combine Stat + Sys limits
  std::vector<TString> sysgraphprefix = {"sysCLgraph_Chi2LD","sysCLgraph_Chi2LD"};
  namecycle = {";2",";2"};
  std::vector<TString> sysnamecycle = {";1",";1"};
  std::array<Double_t,4> leglimits_s = {0.368421,0.140625,0.870927,0.352431};

  CombineCLstatsys(infile,namecycle,sysnamecycle,cat,histtype,histname,histprefix,graphprefix,sysgraphprefix,leglimits_s,lumi,legtitle,outfolder);


}
