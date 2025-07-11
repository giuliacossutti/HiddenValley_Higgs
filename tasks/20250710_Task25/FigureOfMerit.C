/*
Author: Giulia Cossutti

Macro to calculate best signal region like significance and figure of merit of cuts to separate signal from background.
Bkgs are scaled with pythia's calculated cross-sections,
signal is scaled with area = sum of all bkgs areas.
Histograms have been saved in files in 20250705_Task24.

From inside the /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes directory run with
root -l /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250710_Task25/FigureOfMerit.C 

Pre-selection on events:
- At least 2 small jets
- 2 charged leptons with same flavour, opposite charge, | mll - mZ | < tolerance = 10 GeV, ZpT > 20 GeV

Categories of events:
- A: Events with 1 Large Jet
- B: Events with >=2 Large Jets

Selection:
- number of Large Jets >= 1 (to fit in cat. A or B)
- number of b-tagged small jets >=1 (for both A and B events)
- A: DeltaR between first 2 small jets < 2.2
- B: DeltaR between first 2 small jets < 4.8
- A: pT of sum of 2 leptons by Z decay > 96 GeV/c
- A: DeltaR between first 2 small jets <= -0.09 * mjj of first 2 small jets + 8.5
- B: ECF2/pT of Leading Large Jet > 2 GeV/c
*/

//------------------------------------------------------------------------------

#include <numeric>
#include <iostream>

//------------ FUNCTIONS ---------------------
//===================================
// Function to calculate figure of merit with 1D cut
void FOM1D(TFile* infile,const char *histname,Int_t samplenumber,Double_t firstcut,Double_t lastcut,Int_t ncuts,Double_t endcuts, TString outdir){

  if(ncuts <= 1){
    cout << "Use 2 or more cuts." << endl;
    exit(EXIT_FAILURE);
  }

  // Take histograms from files
  std::vector<TH1F*> hist;

  TString histname_i;
  for(Int_t i = 0; i < samplenumber; ++i){
    histname_i.Form("%d",i);
    histname_i = histname + histname_i;
    TH1F *h = (TH1F*)infile->Get(histname_i)->Clone( (histname_i) );
    hist.push_back(h);
  }

  // Calculate figure of merit
  Double_t cut;
  Double_t signal;
  Double_t bkg;
  Double_t bincenter;
  Double_t significance;

  // Canva to draw histograms
  TString canva;
  canva.Form("c%s",histname);

  TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045);
  gPad->SetLeftMargin(0.15);
  //gPad->SetLogy();

  // Graph of figure of merit
  TGraph *fom = new TGraph();
  fom->SetMarkerStyle(kFullCircle);
  fom->SetMarkerColor(kBlack);
  fom->SetLineColor(kRed);
  fom->GetXaxis()->SetTitleSize(.045);
  fom->GetYaxis()->SetTitleSize(.045);
  fom->GetXaxis()->SetTitle("cut");
  fom->GetYaxis()->SetTitle("#frac{S}{#sqrt{Bkg}}");
  TString fomname;
  fomname.Form("FOM_of_%s",histname);
  TString fomtitle;
  fomtitle.Form("FOM of %s, events between cut and %.2f",histname, endcuts);
  fom->SetNameTitle(fomname,fomtitle);

  TH1F *histsignal = new TH1F(fomname + "_sig",fomtitle,ncuts,firstcut-(lastcut-firstcut)/(2.*(ncuts-1)),lastcut+(lastcut-firstcut)/(2.*(ncuts-1)) );
  TH1F *histbkg = new TH1F(fomname + "_bkg",fomtitle,ncuts,firstcut-(lastcut-firstcut)/(2.*(ncuts-1)),lastcut+(lastcut-firstcut)/(2.*(ncuts-1)) );
  histsignal->SetLineColor(kBlue);
  histbkg->SetLineColor(kRed);
  histsignal->GetXaxis()->SetTitleSize(.045);
  histsignal->GetYaxis()->SetTitleSize(.045);
  histsignal->GetXaxis()->SetTitle("cut");
  histsignal->GetYaxis()->SetTitle("Signal or Bkg left (%)");

  // Output file
  TString outfilename;
  outfilename.Form(outdir + "FOM_%s.txt",histname);
  ofstream outfile;
  outfile.open(outfilename);
  outfile << "This file contains significance S/sqrt(Bkg) for the histogram " << histname << endl;
  outfile << "Cuts are between cut and " << endcuts << " with " << ncuts << " cuts between " << firstcut << " and " << lastcut << " included." << endl;
  outfile << "" << endl;

  // Total number of events
  std::vector<Double_t> integrals;
  for(Int_t j = 0; j < hist.size(); ++j){
    integrals.push_back(hist.at(j)->Integral());
  }

  // Loop over cuts
  for(Int_t i = 0; i < ncuts; ++i){
  
    signal = 0.;
    bkg = 0.;

    // Calculate integral cut
    cut = firstcut + i*(lastcut-firstcut)/(ncuts-1);

    // Loop over histograms
    for(Int_t j = 0; j < hist.size(); ++j){
      // Loop over bins of a histogram
      for(Int_t k = 1; k <= hist.at(j)->GetXaxis()->GetNbins(); ++k){
        bincenter = hist.at(j)->GetXaxis()->GetBinCenter(k);
        if( (bincenter < cut && bincenter > endcuts) || (bincenter > cut && bincenter < endcuts) ){
          if(j==0){
            signal += hist.at(j)->GetBinContent(k);
          }else{
            bkg += hist.at(j)->GetBinContent(k);
          }            
        }
      }// end loop over bins of a histogram
    }// end loop over histograms

    // Calculate significance and put it in a graph and in a file with the percentage of considered signal and bkgs
    if(bkg != 0. ){
      significance = signal/sqrt(bkg);
    }else{
      significance = 0.;
    }

    fom->AddPoint(cut,significance);
    outfile <<"Cut:\t"<<cut<<"\tsignificance: "<<significance<<" ,\tLeft "<<100.*signal/integrals.at(0)<<" % signal, \t"<<100.*bkg/accumulate(integrals.begin()+1,integrals.end(),0.)<<" % bkg."<<endl;    

    histsignal->Fill(cut,100.*signal/integrals.at(0) );
    histbkg->Fill(cut,100.*bkg/accumulate(integrals.begin()+1,integrals.end(),0.) );

  }// end loop over cuts

  fom->Draw("APL");
  outfile.close();

  // Save Figure Of Merit in a pdf
  TString outpdf;
  outpdf.Form(outdir + "fom_%s",histname);
  c->SaveAs(outpdf + ".pdf");
  //c->Close();

  // Save histos of signal and bkg left
  TCanvas *csig = new TCanvas(canva + "_sig", canva + "_sig", 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045);
  gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.20);
  //gPad->SetLogy();
  histsignal->SetMaximum(105.);
  histsignal->SetMinimum(0.);

  histsignal->Draw("HIST");
  histbkg->Draw("HIST, same");

  TLegend* legend = new TLegend(0.450501,0.762153,0.549499,0.866319);
  legend->AddEntry(histsignal,"signal");
  legend->AddEntry(histbkg,"background");
  legend->Draw();

  csig->SaveAs(outpdf + "_signal_bkg_left" + ".pdf");
}
//===================================
// Function to calculate figure of merit with 2D rectangular cut
void FOM2D_rectangle(TFile* infile,const char *histname,Int_t samplenumber,Double_t firstcutx,Double_t lastcutx,Int_t ncutsx,Double_t firstcuty,Double_t lastcuty,Int_t ncutsy,TString outdir){

  if(ncutsx <= 1 || ncutsy <= 1){
    cout << "Use 2 or more cuts for each dimension" << endl;
    exit(EXIT_FAILURE);
  }

  // Take histograms from files
  std::vector<TH2F*> hist;

  TString histname_i;
  for(Int_t i = 0; i < samplenumber; ++i){
    histname_i.Form("%d",i);
    histname_i = histname + histname_i;
    TH2F *h = (TH2F*)infile->Get(histname_i)->Clone( (histname_i) );
    hist.push_back(h);
  }

  // Calculate figure of merit
  Double_t cutx;
  Double_t cuty;
  Double_t signal;
  Double_t bkg;
  Double_t bincenterx;
  Double_t bincentery;
  Double_t significance;

  // Canva to draw histograms
  TString canva;
  canva.Form("crectangle_%s",histname);

  TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.20);
  //gPad->SetLogz();

  // 2D histogram of figure of merit
  TString fomname;
  fomname.Form("Rectangular_FOM_of_%s",histname);
  TString fomtitle;
  fomtitle.Form("Rectangular FOM of %s",histname);

  char* axes = strdup(histname);
  char* axesnames1 = strtok(axes,"_");
  char* axesnames2 = strtok(NULL,"_");
  char* axesnames3 = strtok(NULL,"_");
  char* axesnames4 = strtok(NULL,"_");
  char axestitle1[100];
  char axestitle2[100];
  const char* cut = "cut on ";
  strcpy(axestitle1,cut);
  strcat(axestitle1,axesnames2);
  strcpy(axestitle2,cut);
  strcat(axestitle2,axesnames4);
  const char* XTitle = axestitle1;
  const char* YTitle = axestitle2;

  TH2F *fom = new TH2F(fomname,fomtitle,ncutsx,firstcutx-(lastcutx-firstcutx)/(2.*(ncutsx-1)),lastcutx+(lastcutx-firstcutx)/(2.*(ncutsx-1)),ncutsy,firstcuty-(lastcuty-firstcuty)/(2.*(ncutsy-1)),lastcuty+(lastcuty-firstcuty)/(2.*(ncutsy-1)));
  fom->GetXaxis()->SetTitleSize(.045);
  fom->GetYaxis()->SetTitleSize(.045);
  fom->GetZaxis()->SetTitleSize(.045);
  fom->GetXaxis()->SetTitle(XTitle);
  fom->GetYaxis()->SetTitle(YTitle);

  // 2D histograms of signal and bkg left
  TH2F* histsignal = (TH2F*)fom->Clone("Signal_Left_" + fomname);
  TH2F* histbkg = (TH2F*)fom->Clone("Bkg_Left_" + fomname);

  fom->GetZaxis()->SetTitle("#frac{S}{#sqrt{Bkg}}");
  histsignal->GetZaxis()->SetTitle("Signal left (%)");
  histbkg->GetZaxis()->SetTitle("Bkg left (%)");

  // Output file
  TString outfilename;
  outfilename.Form(outdir + "Rectangular_FOM_%s.txt",histname);
  ofstream outfile;
  outfile.open(outfilename);
  outfile << "This file contains significance S/sqrt(Bkg) for the histogram " << histname << endl;
  outfile << "Cuts are between cut and -infinity, with " << ncutsx << " xcuts between " << firstcutx << " and " << lastcutx << " included, " << ncutsy << " ycuts between " << firstcuty << " and " << lastcuty << " included. " <<endl;
  outfile << "XCut is " << XTitle << " , YCut is " << YTitle << endl;
  outfile << "" << endl;

  // Total number of events
  std::vector<Double_t> integrals;
  for(Int_t j = 0; j < hist.size(); ++j){
    integrals.push_back(hist.at(j)->Integral());
  }

  // Loop over cuts
  for(Int_t ix = 0; ix < ncutsx; ++ix){
  for(Int_t iy = 0; iy < ncutsy; ++iy){
  
    signal = 0.;
    bkg = 0.;

    // Calculate cuts
    cutx = firstcutx + ix*(lastcutx-firstcutx)/(ncutsx-1);
    cuty = firstcuty + iy*(lastcuty-firstcuty)/(ncutsy-1);

    // Loop over histograms
    for(Int_t j = 0; j < hist.size(); ++j){
      // Loop over bins of a histogram
      for(Int_t kx = 1; kx <= hist.at(j)->GetXaxis()->GetNbins(); ++kx){
      for(Int_t ky = 1; ky <= hist.at(j)->GetYaxis()->GetNbins(); ++ky){
        bincenterx = hist.at(j)->GetXaxis()->GetBinCenter(kx);
        bincentery = hist.at(j)->GetYaxis()->GetBinCenter(ky);
        if( (bincenterx < cutx) || (bincentery < cuty) ){
          if(j==0){
            signal += hist.at(j)->GetBinContent(kx,ky);
          }else{
            bkg += hist.at(j)->GetBinContent(kx,ky);
          }            
        }
      }// end loop over xbins of a histogram
      }// end loop over ybins of a histogram
    }// end loop over histograms

    // Calculate significance and put it in a histogram and in a file with the percentage of considered signal and bkgs
    if(bkg != 0. ){
      significance = signal/sqrt(bkg);
    }else{
      significance = 0.;
    }

    fom->Fill(cutx,cuty,significance);
    outfile <<"XCut:\t"<<cutx<<" YCut:\t"<<cuty<<"\tsignificance: "<<significance<<" ,\tLeft "<<100.*signal/integrals.at(0)<<" % signal, \t"<<100.*bkg/accumulate(integrals.begin()+1,integrals.end(),0.)<<" % bkg."<<endl;    

    histsignal->Fill(cutx,cuty,100.*signal/integrals.at(0) );
    histbkg->Fill(cutx,cuty,100.*bkg/accumulate(integrals.begin()+1,integrals.end(),0.) );

  }// end loop over ycuts
  }// end loop over xcuts

  fom->Draw("COLZ");
  outfile.close();

  // Save Figure Of Merit in a pdf
  TString outpdf;
  outpdf.Form(outdir + "Rectangular_fom_%s",histname);
  c->SaveAs(outpdf + ".pdf");
  //c->Close();

  // Save histos of signal and bkg left
  TCanvas *csig = new TCanvas(canva + "_sig", canva + "_sig", 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.20);
  //gPad->SetLogz();

  histsignal->Draw("COLZ");
  csig->SaveAs(outpdf + "_signal_left" + ".pdf");

  TCanvas *cbkg = new TCanvas(canva + "_bkg", canva + "_bkg", 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.20);
  //gPad->SetLogz();

  histbkg->Draw("COLZ");
  cbkg->SaveAs(outpdf + "_bkg_left" + ".pdf");
}
//===================================
// Function to calculate figure of merit with 2D line cut
void FOM2D_line(TFile* infile,const char *histname,Int_t samplenumber,Double_t firstcutm,Double_t lastcutm,Int_t ncutsm,Double_t firstcutq,Double_t lastcutq,Int_t ncutsq,TString outdir){

  if(ncutsm <= 1 || ncutsq <= 1){
    cout << "Use 2 or more cuts for each dimension" << endl;
    exit(EXIT_FAILURE);
  }

  // Take histograms from files
  std::vector<TH2F*> hist;

  TString histname_i;
  for(Int_t i = 0; i < samplenumber; ++i){
    histname_i.Form("%d",i);
    histname_i = histname + histname_i;
    TH2F *h = (TH2F*)infile->Get(histname_i)->Clone( (histname_i) );
    hist.push_back(h);
  }

  // Calculate figure of merit
  Double_t cutm;
  Double_t cutq;
  Double_t signal;
  Double_t bkg;
  Double_t bincenterx;
  Double_t bincentery;
  Double_t significance;

  // Canva to draw histograms
  TString canva;
  canva.Form("cline_%s",histname);

  TCanvas *c = new TCanvas(canva, canva, 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.20);
  //gPad->SetLogz();

  // 2D histogram of figure of merit
  char* axes = strdup(histname);
  char* axesnames1 = strtok(axes,"_");
  char* axesnames2 = strtok(NULL,"_");
  char* axesnames3 = strtok(NULL,"_");
  char* axesnames4 = strtok(NULL,"_");

  TString fomname;
  fomname.Form("2D_Linear_FOM_of_%s",histname);
  TString fomtitle;
  fomtitle.Form("2D Linear FOM: Cat %s, %s <= m %s + q",axesnames1,axesnames4,axesnames2);

  const char* XTitle = "m";
  const char* YTitle = "q";

  TH2F *fom = new TH2F(fomname,fomtitle,ncutsm,firstcutm-(lastcutm-firstcutm)/(2.*(ncutsm-1)),lastcutm+(lastcutm-firstcutm)/(2.*(ncutsm-1)),ncutsq,firstcutq-(lastcutq-firstcutq)/(2.*(ncutsq-1)),lastcutq+(lastcutq-firstcutq)/(2.*(ncutsq-1)));
  fom->GetXaxis()->SetTitleSize(.045);
  fom->GetYaxis()->SetTitleSize(.045);
  fom->GetZaxis()->SetTitleSize(.045);
  fom->GetXaxis()->SetTitle(XTitle);
  fom->GetYaxis()->SetTitle(YTitle);

  // 2D histograms of signal and bkg left
  TH2F* histsignal = (TH2F*)fom->Clone("Signal_Left_" + fomname);
  TH2F* histbkg = (TH2F*)fom->Clone("Bkg_Left_" + fomname);

  fom->GetZaxis()->SetTitle("#frac{S}{#sqrt{Bkg}}");
  histsignal->GetZaxis()->SetTitle("Signal left (%)");
  histbkg->GetZaxis()->SetTitle("Bkg left (%)");

  // Output file
  TString outfilename;
  outfilename.Form(outdir + "2D_Linear_FOM_%s.txt",histname);
  ofstream outfile;
  outfile.open(outfilename);
  outfile << "This file contains significance S/sqrt(Bkg) for the histogram " << histname << endl;
  outfile << "Cuts are 2D linear cuts: y <= m x + q  with "<< ncutsm <<" m values between "<< firstcutm <<" and "<< lastcutm <<" included, "<< ncutsq <<" q values between "<< firstcutq <<" and "<< lastcutq <<" included. "<<endl;
  outfile << "x is " << axesnames1 << " , y is " << axesnames3 << endl;
  outfile << "" << endl;

  // Total number of events
  std::vector<Double_t> integrals;
  for(Int_t j = 0; j < hist.size(); ++j){
    integrals.push_back(hist.at(j)->Integral());
  }

  // Loop over cuts
  for(Int_t im = 0; im < ncutsm; ++im){
  for(Int_t iq = 0; iq < ncutsq; ++iq){
  
    signal = 0.;
    bkg = 0.;

    // Calculate m and q values for cuts
    cutm = firstcutm + im*(lastcutm-firstcutm)/(ncutsm-1);
    cutq = firstcutq + iq*(lastcutq-firstcutq)/(ncutsq-1);

    // Loop over histograms
    for(Int_t j = 0; j < hist.size(); ++j){
      // Loop over bins of a histogram
      for(Int_t kx = 1; kx <= hist.at(j)->GetXaxis()->GetNbins(); ++kx){
      for(Int_t ky = 1; ky <= hist.at(j)->GetYaxis()->GetNbins(); ++ky){
        bincenterx = hist.at(j)->GetXaxis()->GetBinCenter(kx);
        bincentery = hist.at(j)->GetYaxis()->GetBinCenter(ky);
        if( bincentery <= bincenterx * cutm + cutq ){
          if(j==0){
            signal += hist.at(j)->GetBinContent(kx,ky);
          }else{
            bkg += hist.at(j)->GetBinContent(kx,ky);
          }            
        }
      }// end loop over xbins of a histogram
      }// end loop over ybins of a histogram
    }// end loop over histograms

    // Calculate significance and put it in a histogram and in a file with the percentage of considered signal and bkgs
    if(bkg != 0. ){
      significance = signal/sqrt(bkg);
    }else{
      significance = 1e8;
    }

    fom->Fill(cutm,cutq,significance);
    outfile <<"mcut:\t"<<cutm<<" qcut:\t"<<cutq<<"\tsignificance: "<<significance<<" ,\tLeft "<<100.*signal/integrals.at(0)<<" % signal, \t"<<100.*bkg/accumulate(integrals.begin()+1,integrals.end(),0.)<<" % bkg."<<endl;    

    histsignal->Fill(cutm,cutq,100.*signal/integrals.at(0) );
    histbkg->Fill(cutm,cutq,100.*bkg/accumulate(integrals.begin()+1,integrals.end(),0.) );

  }// end loop over ycuts
  }// end loop over xcuts

  fom->Draw("COLZ");
  outfile.close();

  // Save Figure Of Merit in a pdf
  TString outpdf;
  outpdf.Form(outdir + "2D_Linear_fom_%s",histname);
  c->SaveAs(outpdf + ".pdf");
  //c->Close();

  // Save histos of signal and bkg left
  TCanvas *csig = new TCanvas(canva + "_sig", canva + "_sig", 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.20);
  //gPad->SetLogz();

  histsignal->Draw("COLZ");
  csig->SaveAs(outpdf + "_signal_left" + ".pdf");

  TCanvas *cbkg = new TCanvas(canva + "_bkg", canva + "_bkg", 0, 0, 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.20);
  //gPad->SetLogz();

  histbkg->Draw("COLZ");
  cbkg->SaveAs(outpdf + "_bkg_left" + ".pdf");

}
//===============================
//---------- END FUNCTIONS -------------------

void FigureOfMerit()
{

  // Open histograms' file
  TString filename = "/gfsvol01/atlas/giuliac/plots_and_outputs/20250710_Task25/area_MergedHistos.root";

  TFile *infile = new TFile(filename);

  // Output directory
  TString outdir = "/gfsvol01/atlas/giuliac/plots_and_outputs/20250710_Task25/FinalSignalRegion_";

  // Number of signal or background samples
  Int_t Nsamples = 5;

  // Calculate 1D FOMs
/*
  FOM1D(infile,"A_DeltaR_between_leading_and_subleading_small_jets",Nsamples,0.4,2.0,9,0.,outdir);
  FOM1D(infile,"B_DeltaR_between_leading_and_subleading_small_jets",Nsamples,0.4,4.6,22,0.,outdir);

  FOM1D(infile,"B_DeltaR_between_leading_and_subleading_large_jets",Nsamples,1.2,3.0,19,0.,outdir+"lowedge_");
  FOM1D(infile,"B_DeltaR_between_leading_and_subleading_large_jets",Nsamples,3.0,6.0,31,10.,outdir+"highedge_");

  FOM1D(infile,"A_2point_Energy_Correlation_div_pT_of_leading_large_jet",Nsamples,0.,40.,21,150.,outdir);
  FOM1D(infile,"B_2point_Energy_Correlation_div_pT_of_leading_large_jet",Nsamples,0.,50.,26,150.,outdir);

  FOM1D(infile,"A_2point_Energy_Correlation_div_pT_of_leading_small_jet",Nsamples,0.,30.,16,150.,outdir);
  FOM1D(infile,"B_2point_Energy_Correlation_div_pT_of_leading_small_jet",Nsamples,0.,30.,16,150.,outdir);

  FOM1D(infile,"A_Invariant_mass_of_all_btagged_small_jets",Nsamples,0.,120.,31,0.,outdir);
  FOM1D(infile,"B_Invariant_mass_of_all_btagged_small_jets",Nsamples,0.,120.,31,0.,outdir);

  FOM1D(infile,"A_Invariant_mass_of_leading_and_subleading_small_jets",Nsamples,0.,120.,31,0.,outdir);
  FOM1D(infile,"B_Invariant_mass_of_leading_and_subleading_small_jets",Nsamples,0.,120.,31,0.,outdir);

  FOM1D(infile,"B_pt_of_first_two_leptons",Nsamples,40.0,200.0,21,600.,outdir);

  FOM1D(infile,"A_DeltaR_between_first_two_leptons",Nsamples,1.0,3.0,21,0.,outdir+"lowedge_");
  FOM1D(infile,"A_DeltaR_between_first_two_leptons",Nsamples,4.0,5.4,15,10.,outdir+"highedge_");
  FOM1D(infile,"B_DeltaR_between_first_two_leptons",Nsamples,1.0,3.0,21,0.,outdir+"lowedge_");
  FOM1D(infile,"B_DeltaR_between_first_two_leptons",Nsamples,4.0,6.0,21,10.,outdir+"highedge_");

  FOM1D(infile,"A_number_of_constituents_of_leading_large_jet",Nsamples,1.5,4.5,4,14.,outdir);
  FOM1D(infile,"B_number_of_constituents_of_leading_large_jet",Nsamples,1.5,4.5,4,14.,outdir);
  FOM1D(infile,"B_number_of_constituents_of_subleading_large_jet",Nsamples,1.5,4.5,4,14.,outdir);

  FOM1D(infile,"A_number_of_small_jets",Nsamples,2.5,5.5,4,14.,outdir);
  FOM1D(infile,"B_number_of_small_jets",Nsamples,2.5,9.5,8,14.,outdir);

  FOM1D(infile,"A_number_of_btagged_small_jets",Nsamples,1.5,4.5,4,14.,outdir);
  FOM1D(infile,"B_number_of_btagged_small_jets",Nsamples,1.5,4.5,4,14.,outdir);

  FOM1D(infile,"A_Prompt_Track_Fraction_of_leading_btagged_small_jet",Nsamples,0.,0.7,36,0.,outdir);
  FOM1D(infile,"B_Prompt_Track_Fraction_of_leading_btagged_small_jet",Nsamples,0.,0.7,36,0.,outdir);

  FOM1D(infile,"A_Prompt_Track_Fraction_of_leading_small_jet",Nsamples,0.,0.7,36,0.,outdir);
  FOM1D(infile,"B_Prompt_Track_Fraction_of_leading_small_jet",Nsamples,0.,0.7,36,0.,outdir);

  FOM1D(infile,"A_Prompt_Track_Fraction_of_subleading_small_jet",Nsamples,0.,0.7,36,0.,outdir);
  FOM1D(infile,"B_Prompt_Track_Fraction_of_subleading_small_jet",Nsamples,0.,0.7,36,0.,outdir);

  FOM1D(infile,"A_Prompt_Track_Fraction_of_leading_large_jet",Nsamples,0.,0.7,36,0.,outdir);
  FOM1D(infile,"B_Prompt_Track_Fraction_of_leading_large_jet",Nsamples,0.,0.7,36,0.,outdir);

  FOM1D(infile,"B_Prompt_Track_Fraction_of_subleading_large_jet",Nsamples,0.,0.7,36,0.,outdir);

  // Calculate 2D Rectangular FOMs
  FOM2D_rectangle(infile,"A_PTF1small_vs_PTF1large",Nsamples,0.,0.3,16,0.,0.3,16, outdir);
  FOM2D_rectangle(infile,"B_PTF1small_vs_PTF1large",Nsamples,0.,0.3,16,0.,0.3,16, outdir);

  FOM2D_rectangle(infile,"A_PTF1bsmall_vs_PTF1large",Nsamples,0.,0.3,16,0.,0.3,16, outdir);
  FOM2D_rectangle(infile,"B_PTF1bsmall_vs_PTF1large",Nsamples,0.,0.3,16,0.,0.3,16, outdir);

  FOM2D_rectangle(infile,"A_PTF1small_vs_PTF2small",Nsamples,0.,0.3,16,0.,0.3,16, outdir);
  FOM2D_rectangle(infile,"B_PTF1small_vs_PTF2small",Nsamples,0.,0.3,16,0.,0.3,16, outdir);
*/
  FOM2D_rectangle(infile,"Chi2HDA_PTF1small_vs_PTF2small",Nsamples,0.05,0.4,8,0.05,0.4,8, outdir);
  FOM2D_rectangle(infile,"Chi2HDB_PTF1small_vs_PTF2small",Nsamples,0.05,0.4,8,0.05,0.4,8, outdir);
/*
  // Calculate 2D Linear FOMs
  FOM2D_line(infile,"A_IMsmall_vs_DeltaRsmall",Nsamples,-0.1,-0.01,10,2.,10.,17, outdir);
  FOM2D_line(infile,"B_IMsmall_vs_DeltaRsmall",Nsamples,-0.1,-0.01,10,3.,10.,15, outdir);
*/
}
