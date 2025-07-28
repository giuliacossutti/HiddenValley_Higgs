/*
Author: Giulia Cossutti

Macro to plot used formulae for systematics evaluation

From inside the /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes directory run with
root -l /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250726_Task28/SysFormulae.C 

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
  header->SetTextSize(0.045);
  gStyle->SetLegendTextSize(0.045);
  legend->SetBorderSize(0);
  legend->SetEntrySeparation(0.01);
}
//=================================
// Function to customize canvas:
void CustomizeCanva(){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleXSize(0.045);
  gStyle->SetTitleYSize(0.045);
  gPad->SetLeftMargin(0.15);
}
//=================================
// b-tagging efficiency formula
Double_t btag(Double_t pt){
  return 0.80*std::tanh(0.003*pt)*(30/(1+0.086*pt));
}
//===================================
// varied b-tagging efficiency formula
Double_t sysbtag(Double_t pt){
  Double_t eff = 0.;
  eff = (pt > 20 && pt <= 30) * ( 0.80*tanh(0.003*pt)*(30/(1+0.086*pt)) ) * ((100.-7.7)/100.) +
        (pt > 30 && pt <= 40) * ( 0.80*tanh(0.003*pt)*(30/(1+0.086*pt)) ) * ((100.-3.0)/100.) +
        (pt > 40 && pt <= 60) * ( 0.80*tanh(0.003*pt)*(30/(1+0.086*pt)) ) * ((100.-1.4)/100.) +
        (pt > 60 && pt <= 85) * ( 0.80*tanh(0.003*pt)*(30/(1+0.086*pt)) ) * ((100.-1.1)/100.) +
        (pt > 85 && pt <= 110) * ( 0.80*tanh(0.003*pt)*(30/(1+0.086*pt)) ) * ((100.-1.0)/100.) +
        (pt > 110 && pt <= 140) * ( 0.80*tanh(0.003*pt)*(30/(1+0.086*pt)) ) * ((100.-1.1)/100.) +
        (pt > 140 && pt <= 175) * ( 0.80*tanh(0.003*pt)*(30/(1+0.086*pt)) ) * ((100.-1.3)/100.) +
        (pt > 175 && pt <= 250) * ( 0.80*tanh(0.003*pt)*(30/(1+0.086*pt)) ) * ((100.-1.5)/100.) +
        (pt > 250)              * ( 0.80*tanh(0.003*pt)*(30/(1+0.086*pt)) ) * ((100.-3.1)/100.);

  return eff;
}
//===================================
// Jet Energy Scale formula
Double_t JES(Double_t pt, Double_t eta){
  return sqrt( pow(3.0 - 0.2*(abs(eta)),2) / pt + 1.0 );
}
//===================================
// varied Jet Energy Scale formula
Double_t sysJES(Double_t pt, Double_t eta){
  return sqrt( pow(3.0 - 0.2*(abs(eta)),2) / pt + 1.0 ) * (96./100.);
}
//===================================
//---------- END FUNCTIONS -------------------


void SysFormulae()
{
  TString outdir = "/gfsvol01/atlas/giuliac/plots_and_outputs/20250726_Task28/";

  // b-tagging efficiency
  TF1* btag_nom = new TF1("btag_nom","btag(x)",20.,600.);
  TF1* btag_sys = new TF1("btag_sys","sysbtag(x)",20.,600.);
  btag_sys->SetNpx(1000);
  btag_nom->SetLineColor(kGreen-2);
  btag_sys->SetLineColor(kBlack);

  // Canva
  TString canva = "cbtag";
  TCanvas *c1 = new TCanvas(canva, canva, 0, 0, 800, 600);
  CustomizeCanva();
  gPad->SetLogx();

  btag_nom->Draw();
  btag_sys->Draw("sames");

  btag_nom->GetYaxis()->SetRangeUser(0.45,0.95);
  btag_nom->GetXaxis()->SetMoreLogLabels();
  btag_nom->GetXaxis()->SetNoExponent();	
  btag_nom->GetYaxis()->SetTitle("#varepsilon_{b}");
  btag_nom->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  // Legend
  std::array<Double_t,4> leglimits = {0.299499,0.616319,0.799499,0.847222};

  TLegend* legend = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  TString head;
  head = TString::Format("b-tagging Efficiency");

  legend->SetHeader(head,"C");
  CustomizeLeg(legend);

  // Add the right entry to the legend
  legend->AddEntry(btag_nom,"Nominal");
  legend->AddEntry(btag_sys,"1#sigma Systematic Deviation");
  legend->Draw();

  c1->SaveAs(outdir + "btageff.pdf");




  // Jet Energy Scale
  TF2* JES_nom = new TF2("JES_nom","JES(x,y)",20.,600.,-2.6,2.6);
  TF2* JES_sys = new TF2("JES_sys","sysJES(x,y)",20.,600.,-2.6,2.6);

  // Projections on pt and eta axes
  Double_t eta = 0.0;
  TF12* JES_nom_pt = new TF12("JES_nom_pt",JES_nom,eta,"x");
  TF12* JES_sys_pt = new TF12("JES_sys_pt",JES_sys,eta,"x");
  JES_nom_pt->SetLineColor(kGreen-2);
  JES_sys_pt->SetLineColor(kBlack);
  JES_sys_pt->SetNpx(1000);

  Double_t pt = 60.;
  TF12* JES_nom_eta = new TF12("JES_nom_eta",JES_nom,pt,"y");
  TF12* JES_sys_eta = new TF12("JES_sys_eta",JES_sys,pt,"y");
  JES_nom_eta->SetLineColor(kGreen-2);
  JES_sys_eta->SetLineColor(kBlack);
  JES_sys_eta->SetNpx(1000);


  // Canva
  canva = "cJESpt";
  TCanvas *c2 = new TCanvas(canva, canva, 0, 0, 800, 600);
  CustomizeCanva();
  gPad->SetLogx();

  JES_nom_pt->Draw();
  JES_sys_pt->Draw("sames");

  JES_nom_pt->GetXaxis()->SetMoreLogLabels();
  JES_nom_pt->GetXaxis()->SetNoExponent();	
  JES_nom_pt->GetYaxis()->SetTitle("JES");
  JES_nom_pt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  JES_nom_pt->GetYaxis()->SetRangeUser(0.95,1.25);

  // Legend
  TLegend* legendjp = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  head = TString::Format("Jet Energy Scale Formula, #eta = %.1f", eta);

  legendjp->SetHeader(head,"C");
  CustomizeLeg(legendjp);

  // Add the right entry to the legend
  legendjp->AddEntry(JES_nom_pt,"Nominal");
  legendjp->AddEntry(JES_sys_pt,"1#sigma Systematic Deviation");
  legendjp->Draw();

  c2->SaveAs(outdir + "JESpt.pdf");

  // Canva
  canva = "cJESeta";
  TCanvas *c3 = new TCanvas(canva, canva, 0, 0, 800, 600);
  CustomizeCanva();
  //gPad->SetLogx();

  JES_nom_eta->Draw();
  JES_sys_eta->Draw("sames");

  //JES_nom_pt->GetXaxis()->SetMoreLogLabels();
  //JES_nom_pt->GetXaxis()->SetNoExponent();	
  JES_nom_eta->GetYaxis()->SetTitle("JES");
  JES_nom_eta->GetXaxis()->SetTitle("#eta");
  JES_nom_eta->GetYaxis()->SetRangeUser(0.95,1.25);

  // Legend
  TLegend* legendje = new TLegend(leglimits[0],leglimits[1],leglimits[2],leglimits[3]);

  head = TString::Format("Jet Energy Scale Formula, p_{T} = %.1f GeV/c", pt);

  legendje->SetHeader(head,"C");
  CustomizeLeg(legendje);

  // Add the right entry to the legend
  legendje->AddEntry(JES_nom_eta,"Nominal");
  legendje->AddEntry(JES_sys_eta,"1#sigma Systematic Deviation");
  legendje->Draw();

  c3->SaveAs(outdir + "JESeta.pdf");




  // Tracking efficiency

  // pt bins
  Int_t nbinsx = 3;
  const Double_t xbins[] = {0.,0.1,1.0,2.0};

  // eta bins
  Int_t nbinsy = 3;
  const Double_t ybins[] = {0.,1.5,2.5,3.5};

  // histograms of efficiencies
  TH2F* etrk_nom = new TH2F("etrk_nom","",nbinsx,xbins,nbinsy,ybins);
  //TH2F* etrk_sys = new TH2F("etrk_sys","",nbinsx,xbins,nbinsy,ybins);
  TH2F* etrk_sys = (TH2F*)etrk_nom->Clone("etrk_sys");

  // Set efficiency values
  for(Int_t ky = 1; ky <= etrk_nom->GetYaxis()->GetNbins(); ++ky){
    etrk_nom->SetBinContent(1,ky,0.);
    etrk_sys->SetBinContent(1,ky,0.);
  }
  for(Int_t kx = 1; kx <= etrk_nom->GetXaxis()->GetNbins(); ++kx){
    etrk_nom->SetBinContent(kx,3,0.);
    etrk_sys->SetBinContent(kx,3,0.);
  }

  etrk_nom->SetBinContent(2,1,0.70);
  etrk_nom->SetBinContent(2,2,0.60);
  etrk_nom->SetBinContent(3,1,0.95);
  etrk_nom->SetBinContent(3,2,0.85);

  etrk_sys->SetBinContent(2,1,0.70 * 0.90);
  etrk_sys->SetBinContent(2,2,0.60 * 0.92);
  etrk_sys->SetBinContent(3,1,0.95 * 0.90);
  etrk_sys->SetBinContent(3,2,0.85 * 0.92);

  Double_t Max = 1.;
  Double_t Min = 0.;

  std::array<Double_t,4> leglimits2 = {0.299499,0.692708,0.799499,0.873264};

  // Canva
  canva = "cetrkn";
  TCanvas *c4 = new TCanvas(canva, canva, 0, 0, 800, 600);
  CustomizeCanva();
  gPad->SetRightMargin(0.15);
  //gPad->SetLogx();
  gStyle->SetPaintTextFormat(".2f");

  etrk_nom->SetMarkerSize(2.);
  etrk_nom->Draw("text COLZ0");
  etrk_nom->SetMaximum(Max);
  etrk_nom->SetMinimum(Min);

  etrk_nom->GetZaxis()->SetTitleSize(.045);
  etrk_nom->GetYaxis()->SetTitle("|#eta|");
  etrk_nom->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  etrk_nom->GetZaxis()->SetTitle("#varepsilon_{trk}");

  // Legend
  TLegend* legendtn = new TLegend(leglimits2[0],leglimits2[1],leglimits2[2],leglimits2[3]);

  head = TString::Format("Charged Hadron Tracking Efficiency");

  legendtn->SetHeader(head,"C");
  CustomizeLeg(legendtn);

  // Add the right entry to the legend
  legendtn->AddEntry(etrk_nom,"Nominal");
  legendtn->Draw();

  c4->SaveAs(outdir + "etrknom.pdf");


  // Canva
  canva = "cetrks";
  TCanvas *c5 = new TCanvas(canva, canva, 0, 0, 800, 600);
  CustomizeCanva();
  gPad->SetRightMargin(0.15);
  //gPad->SetLogx();
  gStyle->SetPaintTextFormat(".2f");
 
  etrk_sys->SetMarkerSize(2.);
  etrk_sys->Draw("text COLZ0");
  etrk_sys->SetMaximum(Max);
  etrk_sys->SetMinimum(Min);

  etrk_sys->GetZaxis()->SetTitleSize(.045);
  etrk_sys->GetYaxis()->SetTitle("|#eta|");
  etrk_sys->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  etrk_sys->GetZaxis()->SetTitle("#varepsilon_{trk}");

  // Legend
  TLegend* legendts = new TLegend(leglimits2[0],leglimits2[1],leglimits2[2],leglimits2[3]);

  head = TString::Format("Charged Hadron Tracking Efficiency");

  legendts->SetHeader(head,"C");
  CustomizeLeg(legendts);

  // Add the right entry to the legend
  legendts->AddEntry(etrk_sys,"1#sigma Systematic Deviation");
  legendts->Draw();

  c5->SaveAs(outdir + "etrksys.pdf");

}
