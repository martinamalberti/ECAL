//****** simple macro to compute PU weights ******
{
  //*** mc file
  TChain *ntu_MC = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
  ntu_MC->Add("/media/DATA/ALCARECO/DYSummer12/DYJets-Summer12.root");
  TH1F *hmc  = new TH1F("hmc","hmc",60,0,60);
  ntu_MC->Draw("PUit_TrueNumInteractions >> hmc","","GOFF");
 
  //*** data file 
  TFile *fda = TFile::Open("/media/DATA/ALCARECO/MyDataPileupHistogram_minBiasXsec69400_190456-195016.root");
  TH1F *hdata = (TH1F*)fda->Get("pileup");
 
  //*** compute weights
  TH1F *hweights = (TH1F*)hdata->Clone("hweights");
  hweights->Divide(hdata,hmc,1./hdata->GetSumOfWeights(),1./hmc->GetSumOfWeights());

  TFile *fout = new TFile("./PUweights_2012_DYJetsToLL_Summer12_S9_minBiasXsec69400_190456-195016.root","create");
  hweights->Write("hweights");



}
