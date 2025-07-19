/*
Author: Giulia Cossutti

Macro to calculate 1 sigma impact of systematics to background samples.
Histograms of bin-per-bin deviations in the histogram chosen for 95% CL limits calculation are saved in files.

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
root -l /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250718_Task27/SysSigma.C 

*/

//------------------------------------------------------------------------------

//------------ FUNCTIONS ---------------------
//===================================
// Function to calculate difference between two histograms from different files and write it in another file
void Diff1Sigma(std::vector<TFile*> fnominal, std::vector<TFile*> fsys, std::vector<TString> bkgtag, TString systag, const char *histname, std::vector<Double_t> weights){

  // Take histograms from files
  std::vector<TH2F*> histnom;
  std::vector<TH2F*> histsys;

  for(Int_t i = 0; i < bkgtag.size(); ++i){
    TH2F *hn = (TH2F*)fnominal.at(i)->Get(histname)->Clone( TString::Format("nom_%s_%d", histname, i).Data() );
    histnom.push_back(hn);
    TH2F *hs = (TH2F*)fsys.at(i)->Get(histname)->Clone( TString::Format("sys_%s_%d" + systag , histname , i).Data() );
    histsys.push_back(hs);
  }

  // Scale histograms by Integrated Luminosity, Cross-section and number of generated events
  for(Int_t i = 0; i < bkgtag.size(); ++i){
    histnom.at(i)->Sumw2();
    histnom.at(i)->Scale(weights.at(i));
    histsys.at(i)->Sumw2();
    histsys.at(i)->Scale(weights.at(i));

    // Write on file
    histnom.at(i)->Write("",TObject::kOverwrite);
    histsys.at(i)->Write("",TObject::kOverwrite);
  }

  // Histogram of sum of bkgs
  TH2F *histnom_bkg = (TH2F*)histnom.at(0)->Clone( TString::Format("bkg_nom_%s", histname).Data() );
  TH2F *histsys_bkg = (TH2F*)histsys.at(0)->Clone( TString::Format("bkg_sys_%s" + systag , histname).Data() );
  for(Int_t i = 0; i < bkgtag.size(); ++i){
    histnom_bkg->Add(histnom.at(i));
    histsys_bkg->Add(histsys.at(i));
  }

  // Difference of nominal and sys histograms
  TH2F *histdiff = (TH2F*)histnom_bkg->Clone( TString::Format( "%s_sys" + systag, histname).Data() );
  histdiff->Add(histsys_bkg, -1.);

  // Write on file
  histnom_bkg->Write("",TObject::kOverwrite);
  histsys_bkg->Write("",TObject::kOverwrite);
  histdiff->Write("",TObject::kOverwrite);
}
//===================================
// Function to calculate difference between two histograms and write it in another file
void Diff1Sigma_weights(std::vector<TFile*> fnominal, std::vector<TString> bkgtag, TString systag, const char *histname, std::vector<Double_t> weights, std::vector<Double_t> sysweights){

  // Take histograms from files
  std::vector<TH2F*> histnom;
  std::vector<TH2F*> histsys;

  for(Int_t i = 0; i < bkgtag.size(); ++i){
    TH2F *hn = (TH2F*)fnominal.at(i)->Get(histname)->Clone( TString::Format("nom_%s_%d" , histname, i).Data() );
    histnom.push_back(hn);
    TH2F *hs = (TH2F*)fnominal.at(i)->Get(histname)->Clone( TString::Format("sys_%s_%d" + systag , histname , i).Data() );
    histsys.push_back(hs);
  }

  // Scale histograms by Integrated Luminosity, Cross-section and number of generated events
  for(Int_t i = 0; i < bkgtag.size(); ++i){
    histnom.at(i)->Sumw2();
    histnom.at(i)->Scale(weights.at(i));
    histsys.at(i)->Sumw2();
    histsys.at(i)->Scale(sysweights.at(i));

    // Write on file
    histnom.at(i)->Write("",TObject::kOverwrite);
    histsys.at(i)->Write("",TObject::kOverwrite);
  }

  // Histogram of sum of bkgs
  TH2F *histnom_bkg = (TH2F*)histnom.at(0)->Clone( TString::Format("bkg_nom_%s", histname).Data() );
  TH2F *histsys_bkg = (TH2F*)histsys.at(0)->Clone( TString::Format("bkg_sys_%s" + systag , histname).Data() );
  for(Int_t i = 0; i < bkgtag.size(); ++i){
    histnom_bkg->Add(histnom.at(i));
    histsys_bkg->Add(histsys.at(i));
  }

  // Difference of nominal and sys histograms
  TH2F *histdiff = (TH2F*)histnom_bkg->Clone( TString::Format( "%s_sys" + systag, histname).Data() );
  histdiff->Add(histsys_bkg, -1.);

  // Write on file
  histnom_bkg->Write("",TObject::kOverwrite);
  histsys_bkg->Write("",TObject::kOverwrite);
  histdiff->Write("",TObject::kOverwrite);
}
//===================================
//---------- END FUNCTIONS -------------------


void SysSigma()
{

  gSystem->Load("libDelphes");
  
  // Import histograms 
  TString dir = "/gfsvol01/atlas/giuliac/plots_and_outputs/20250718_Task27/";

  // File names and tags
  std::vector<TString> bkgtag = {"bkg_HZ_SM","bkg_ZZ","bkg_Zjets_g","bkg_Zjets_q"};
  std::vector<TString> systag = {"_JES", "_etrk", "_eb"};

  std::vector<TFile*> fnominal;
  std::vector<std::vector<TFile*>> fsys;

  // Nominal files
  for(Int_t i = 0; i < bkgtag.size(); ++i){
    TFile *fn = new TFile(dir + bkgtag.at(i) + ".root");
    fnominal.push_back(fn);
  }

  // Files with systematic number j applied
  for(Int_t j = 0; j < systag.size(); ++j){
    std::vector<TFile*> fsj;
    for(Int_t i = 0; i < bkgtag.size(); ++i){
      TFile *fs = new TFile(dir + bkgtag.at(i) + systag.at(j) + ".root");
      fsj.push_back(fs);
    }    
    fsys.push_back(fsj);
  }

  // Control size of inputs
  if( !( (bkgtag.size() == fnominal.size()) && (systag.size() == fsys.size()) )  ){
    cout << "Error in SysSigma: number of files not defined" << endl;
    cout << "fnominal size: " << fnominal.size() << "fsys size: " << fsys.size() << endl;
    exit(EXIT_FAILURE);
  }
  
  // Integrated Luminosity [pb^-1]
  Double_t lumi = 360000;

  // Cross-section [pb] 
  //Double_t xsec_bkgZq = 1.250e4;
  //Double_t xsec_bkgZg = 7.759e3;
  Double_t xsec_bkgZq = 1.671e3;
  Double_t xsec_bkgZg = 4.374e2;
  Double_t xsec_bkgZZ = 7.688e-1;
  Double_t xsec_bkgHZ = 6.455e-2;

  // Number of MC events
  Double_t ntot_bkgZq = 5000000;
  Double_t ntot_bkgZg = 5000000;
  Double_t ntot_bkgZZ = 500000;
  Double_t ntot_bkgHZ = 500000;

  // Calculate weights
  std::vector<Double_t> xsec = {xsec_bkgHZ,xsec_bkgZZ,xsec_bkgZg,xsec_bkgZq};
  std::vector<Double_t> ntot = {ntot_bkgHZ,ntot_bkgZZ,ntot_bkgZg,ntot_bkgZq};

  // Control size of inputs
  if( !( (bkgtag.size() == xsec.size()) && (bkgtag.size() == ntot.size()) )  ){
    cout << "Error in SysSigma: number of samples not defined" << endl;
    cout << "bkgtag size: " << bkgtag.size() << "   xsec size: " << xsec.size() << endl;
    cout << "ntot size: " << ntot.size() << endl;
    exit(EXIT_FAILURE);
  }

  std::vector<Double_t> weights;

  for(Int_t i = 0; i < bkgtag.size(); ++i){
    weights.push_back( lumi*xsec.at(i) / ntot.at(i) );
  }

  // write all scaled histograms on file
  TString outname = dir + "sys.root";
  TFile *f = new TFile(outname,"RECREATE");

  // Event categories
  std::vector<TString> cat = {"A_","B_"};

  // Calculate 1 sigma differences for the systag systematics
  for(Int_t i = 0; i < cat.size(); ++i){
    for(Int_t j = 0; j < systag.size(); ++j){
      Diff1Sigma(fnominal, fsys.at(j), bkgtag, systag.at(j), "Chi2LD" + cat.at(i) + "PTF1small_vs_PTF2small", weights);
    }
  }

  // Calculate 1 sigma differences for the Z+jets xsec systematics
  Double_t sysxsec_bkgZq = xsec_bkgZq * ((100.-5.6)/100.);
  Double_t sysxsec_bkgZg = xsec_bkgZg * ((100.-5.6)/100.);

  std::vector<Double_t> sysxsec = {xsec_bkgHZ,xsec_bkgZZ,sysxsec_bkgZg,sysxsec_bkgZq};

  // Control size of inputs
  if( !( (bkgtag.size() == sysxsec.size()) && (bkgtag.size() == ntot.size()) )  ){
    cout << "Error in SysSigma: number of samples not defined in sys_xsec" << endl;
    cout << "bkgtag size: " << bkgtag.size() << "   sysxsec size: " << sysxsec.size() << endl;
    cout << "ntot size: " << ntot.size() << endl;
    exit(EXIT_FAILURE);
  }

  std::vector<Double_t> sysweights;

  for(Int_t i = 0; i < bkgtag.size(); ++i){
    sysweights.push_back( lumi*sysxsec.at(i) / ntot.at(i) );
  }

  for(Int_t i = 0; i < cat.size(); ++i){
    Diff1Sigma_weights(fnominal, bkgtag, "_xsec", "Chi2LD" + cat.at(i) + "PTF1small_vs_PTF2small", weights, sysweights);
  }

}
