// per compilare: g++ -Wall -o LocalCorrectionEBDoubleCBFit `root-config --cflags --glibs` -lboost_regex -lboost_program_options LocalCorrectionEBDoubleCBFit.cpp ../CommonTools/templateUtils.cc

#include "../CommonTools/histoFunc.h"
#include "../CommonTools/TPileupReweighting.h"
#include "../CommonTools/templateUtils.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>
#include <typeinfo>

#include "boost/program_options.hpp"


#define XTALWIDTH 0.01745329
#define PI        3.1415926536 

using namespace std;
using namespace boost;
namespace po = boost::program_options;

//--- options parsing
string outfilename_;
float etaMax_ = 1.442;
float r9min_  = 0.;
float r9max_  = 999.;  
float bcNmin_ = 0;
float bcNmax_ = 999;

bool useW_         = false;
bool useZ_         = false;
bool usePUweights_ = false;
bool useOddCry_    = false;
bool useEvenCry_   = false;
bool useLocalPhi_  = false;
bool applyLocalCorrections_ = false;

void OptionParser(int argc, char *argv[]){
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                        "Show help")
    ("outfilename,o", po::value<string>(&outfilename_)->default_value("graphs.root"), "Output file name")
    ("minR9", po::value<float>(&r9min_)->default_value(0.),                            "Min R9")
    ("maxR9", po::value<float>(&r9max_)->default_value(999),                           "Max R9")
    ("minNBC", po::value<float>(&bcNmin_)->default_value(0.),                          "Min number of basic clusters")
    ("maxNBC", po::value<float>(&bcNmax_)->default_value(999),                         "Max number of basic clusers")
    ("useW",                                                                           "Use W events")
    ("useZ",                                                                           "Use Z events")
    ("usePUweights",                                                                   "Use pile-up weights in MC")
    ("useOddCry",                                                                      "Use odd crystals only")
    ("useEvenCry",                                                                     "Use even crystals only")
    ("useLocalPhi",                                                                    "Analyze vs Local Phi")
    ("applyLocalCorrections",                                                          "Apply local corrections offline")
    ;                                                                             
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")){ cout << desc << endl; exit(1);}
  if (vm.count("useW"))             useW_=true;
  if (vm.count("useZ"))             useZ_=true;
  if (vm.count("usePUweights"))     usePUweights_=true;
  if (vm.count("useOddCry"))        useOddCry_=true;
  if (vm.count("useEvenCry"))       useEvenCry_=true;
  if (vm.count("useLocalPhi"))      useLocalPhi_=true;
  if (vm.count("applyLocalCorrections")) applyLocalCorrections_=true;
}
//------------------------------------------------------------------

//----crystal in gap -----------------------------------------------
bool IsEtaGap(float eta){
  float feta = fabs(eta);
  if( fabs(feta - 0 )<2) return true;
  if( fabs(feta - 25)<2) return true;
  if( fabs(feta - 45)<2) return true;
  if( fabs(feta - 65)<2) return true;
  if( fabs(feta - 85)<2) return true;
  return false;
}
//------------------------------------------------------------------

//----template index -----------------------------------------------
int templIndex(float eta){
    float feta = fabs(eta);
    if (feta <= 25)            {return 0;}
    if (feta> 25 && feta <= 45){return 1;}
    if (feta> 45 && feta <= 65){return 2;}
    if (feta> 65 && feta <= 85){return 3;}

    return -1;
}
//------------------------------------------------------------------


//---- fill histograms -----------------------------------------------
void FillHistograms(float energy, float energy_regr, float p, TH1F* h_template, TH1F* hc_template, TH1F* h, TH1F* hc, float weight, int mod, int bin, int nBins){

  //-- fill templates for each mod                                                                                                                                                                                      
  float correction = energy_regr/energy;
  h_template  -> Fill(energy/p,weight);
  hc_template -> Fill(energy_regr/p,weight);

  //-- fill MC histos in eta bins
  if (bin>nBins-1 || bin < 0 ) {                                                                                                                                                                                      
    std::cout << "Error in bins: " << bin << " . Max number of bins is"<< nBins<< std::endl;
    return;
  }                                                                                                                                                                                                                   

  h  -> Fill(energy/p,weight);                                                                                                                                                                            
  hc -> Fill(energy_regr/p,weight);   
}
//------------------------------------------------------------------

void ActivateBranches(TTree *tree, std::vector<string> list, bool isMC){
  tree->SetBranchStatus("*",0);
  for (int i = 0; i<list.size(); i++){
    if (!isMC && list[i]=="PUit_NumInteractions") 
      continue;
    tree->SetBranchStatus(list[i].c_str(),1);
  }
}


//------------------------------------------------------------------
void DoFit(TH1F *h, TF1 *fun, float &scale, float &err){

  float initialPars[7]={100.,1.,0.05,0.8,2.,1.5,3.5};

  fun -> SetParameter(0, h->GetMaximum());
  for (int ipar = 1; ipar < fun->GetNumberFreeParameters(); ipar++){
    fun -> SetParameter(ipar,initialPars[ipar]);
  }
  int fitStatus = h -> Fit(fun->GetName(), "SMRQ");
  int nTrials = 0;
  while( (fitStatus != 1) && (nTrials < 5) )
    {
      fun -> SetParameter(4, initialPars[4]+0.4*nTrials );
      fun -> SetParameter(6, initialPars[4]+0.4*nTrials );
      fitStatus = h -> Fit(fun->GetName(), "SMRQL+");
      if( fitStatus == 1 ) break;
      ++nTrials;
    }
  scale = fun->GetParameter(1);
  err   = fun->GetParError(1);
}



//---- MAIN PROGRAM -----------------------------------------------
int main(int argc, char** argv)
{
  TF1 *feta[4];
  feta[0] = new TF1("feta0","[0]+([1]*(x)+[2]*pow(x,2))",-1,1.);
  feta[1] = new TF1("feta1","[0]+([1]*(x)+[2]*pow(x,2))",-1,1.);
  feta[2] = new TF1("feta2","[0]+([1]*(x)+[2]*pow(x,2))",-1,1.);
  feta[3] = new TF1("feta3","[0]+([1]*(x)+[2]*pow(x,2))",-1,1.);
  feta[0]->SetParameters(1.00603 , 0.00300789 , -0.0667232);
  feta[1]->SetParameters(1.00655 , 0.00386189 , -0.073931 );
  feta[2]->SetParameters(1.00634 , 0.00631341 , -0.0764134);
  feta[3]->SetParameters(1.00957 , 0.0113306 ,  -0.123808 );
 
  TF1 *fphi[4] ;
  fphi[0] = new TF1("fphi0","[0]+([1]*(x)+[2]*pow(x,2))",-1,1.);
  fphi[1] = new TF1("fphi1","[0]+([1]*(x)+[2]*pow(x,2))",-1,1.);
  fphi[2] = new TF1("fphi2","[0]+([1]*(x)+[2]*pow(x,2))",-1,1.);
  fphi[3] = new TF1("fphi3","[0]+([1]*(x)+[2]*pow(x,2))",-1,1.);
  fphi[0]->SetParameters(1.00403 , -0.0012733 , -0.042925);
  fphi[1]->SetParameters(1.00394 , -0.00137567 , -0.0416698);
  fphi[2]->SetParameters(1.00298 , -0.00111589 , -0.0320377);
  fphi[3]->SetParameters(1.00269 , -0.00153347 , -0.0296769);


  //-- parse options
  OptionParser(argc,argv);
  std::cout<< "Output file name : " << outfilename_.c_str()<<std::endl;
  std::cout<< "Min R9  : " << r9min_<<std::endl;
  std::cout<< "Max R9  : " << r9max_<<std::endl;
  std::cout<< "Min nBC : " << bcNmin_<<std::endl;
  std::cout<< "Max nBC : " << bcNmax_<<std::endl;
  std::cout<< "Using W events : " << useW_<<std::endl;
  std::cout<< "Using Z events : " << useZ_<<std::endl;
    
  //---- PU weights for MC
  TPileupReweighting *puReweighting;
  if (useW_) puReweighting = new TPileupReweighting("~/public/Pileup/PUweights_DYJetsToLL_Summer12_DR53X-PU_S10_minBiasXsec69400_corr_observed_Run2012ABCD.root","hweights"); 
  if (useZ_) puReweighting = new TPileupReweighting("~/public/Pileup/PUweights_DYJetsToLL_Summer12_DR53X-PU_RD1_minBiasXsec69400_pixelcorr_observed_2012ABCD.root","hweights"); 
    
  //---- NTUPLES
  TChain *ntu_MC = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
  TChain *ntu_Data = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
    
  //---- MC Summer 2012 
  if (useW_){
    //    ntu_MC->Add("~/eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/MC_53X/DYJets-Summer12-START53-noSkim/simpleNtupleEoverP_mc_*.root");
    ntu_MC->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM.root");
  }
  if (useZ_){
    //    ntu_MC->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/MC_53X/DYJets-Summer12-START53-noSkim/simpleNtupleEoverP_mc_*.root");
    ntu_MC->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_RD1_START53_V7N-v1_AODSIM.root");
  }

  //---- DATA
  if (useW_){
    // -- rereco 13Jul
    //ntu_Data->Add("~/eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/ReReco/13Jul2012_v2/SingleElectron_Run2012A-WElectron-13Jul2012-v1_USER/SingleElectron_Run2012A-WElectron-13Jul2012-v1_USER_ts.root");
    //ntu_Data->Add("~/eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/ReReco/13Jul2012_v2/SingleElectron_Run2012B-WElectron-13Jul2012-v1_USER/SingleElectron_Run2012B-WElectron-13Jul2012-v1_USER_ts.root");
    // -- rereco 22 Jan
    ntu_Data->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/ReReco/22Jan2013-v1_6/SingleElectron_Run2012A-22Jan2013-v1_AOD/simpleNtuple_*.root");
    ntu_Data->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/ReReco/22Jan2013-v1_6/SingleElectron_Run2012B-22Jan2013-v1_AOD/simpleNtuple_*.root");
    ntu_Data->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/ReReco/22Jan2013-v1_6/SingleElectron_Run2012C-22Jan2013-v1_AOD/simpleNtuple_*.root");
    ntu_Data->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/ReReco/22Jan2013-v1_6/SingleElectron_Run2012D-22Jan2013-v1_AOD/simpleNtuple_*.root");
  }
  if (useZ_){
    //ntu_Data->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/ReReco/13Jul2012_v2/DoubleElectron_Run2012A-ZElectron-13Jul2012-v1_RAW-RECO/simpleNtuple*.root");
    //ntu_Data->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/ReReco/13Jul2012_v2/DoubleElectron_Run2012B-ZElectron-13Jul2012-v1_RAW-RECO/simpleNtuple*.root");
    ntu_Data->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/ReReco/22Jan2013-v1_6/DoubleElectron_Run2012A-22Jan2013-v1_AOD.root");
    ntu_Data->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/ReReco/22Jan2013-v1_6/DoubleElectron_Run2012B-22Jan2013-v1_AOD.root");
    ntu_Data->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/ReReco/22Jan2013-v1_6/DoubleElectron_Run2012C-22Jan2013-v1_AOD.root");
    ntu_Data->Add("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalMIBI/NTUPLES_EOverP/ReReco/22Jan2013-v1_6/DoubleElectron_Run2012D-22Jan2013-v1_AOD.root");
  }
 
  std::cout << "     MC    : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;
  std::cout << "     Data  : " << ntu_Data->GetEntries() << " entries in  Data  sample" << std::endl;

  //---- Observables
  int isZ;
  int npu;
  float npuTrue;
  float PV_z;
  int PV_n;
  float ele1_EoP, ele1_scEta, ele1_scPhi;
  float ele1_scE3x3, ele1_scERaw, ele1_scE, ele1_scE_regression, ele1_seedE;  
  float ele1_charge, ele1_scLocalEta, ele1_scLocalPhi;
  float ele1_R9;
  float ele1_tkP;
  int ele1_bcN;
  float ele1_eSeedBC;
  float ele2_EoP, ele2_scEta, ele2_scPhi;
  float ele2_scE3x3, ele2_scERaw, ele2_scE, ele2_scE_regression, ele2_seedE;  
  float ele2_charge, ele2_scLocalEta, ele2_scLocalPhi;
  float ele2_R9;
  float ele2_tkP;
  int ele2_bcN;
  float ele2_eSeedBC;

  std::vector<std::string> branchesToRead;
  branchesToRead.push_back("PUit_NumInteractions");
  branchesToRead.push_back("isZ");
  branchesToRead.push_back("ele1_scEta");
  branchesToRead.push_back("ele1_scPhi");
  branchesToRead.push_back("ele1_EOverP");
  branchesToRead.push_back("ele1_tkP");
  branchesToRead.push_back("ele1_e3x3");
  branchesToRead.push_back("ele1_scERaw");
  branchesToRead.push_back("ele1_scE");
  branchesToRead.push_back("ele1_scE_regression");
  branchesToRead.push_back("ele1_eSeedBC");
  branchesToRead.push_back("ele1_charge");
  branchesToRead.push_back("ele1_scLocalPhi"); 
  branchesToRead.push_back("ele1_scLocalEta"); 
  branchesToRead.push_back("ele1_bcN"); 
  if (useZ_){
    branchesToRead.push_back("ele2_scEta");
    branchesToRead.push_back("ele2_scPhi");
    branchesToRead.push_back("ele2_EOverP");
    branchesToRead.push_back("ele2_tkP");
    branchesToRead.push_back("ele2_e3x3");
    branchesToRead.push_back("ele2_scERaw");
    branchesToRead.push_back("ele2_scE");
    branchesToRead.push_back("ele2_scE_regression");
    branchesToRead.push_back("ele2_eSeedBC");
    branchesToRead.push_back("ele2_charge");
    branchesToRead.push_back("ele2_scLocalPhi"); 
    branchesToRead.push_back("ele2_scLocalEta"); 
    branchesToRead.push_back("ele2_bcN"); 
  }


  //---- Set branch addresses for MC  
  ActivateBranches(ntu_MC,branchesToRead,true);
  ntu_MC->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_MC->SetBranchAddress("isZ", &isZ);
  //ntu_MC->SetBranchAddress("PV_z", &PV_z);
  //ntu_MC->SetBranchAddress("PV_n", &PV_n);
  ntu_MC->SetBranchAddress("ele1_scEta", &ele1_scEta);
  ntu_MC->SetBranchAddress("ele1_scPhi", &ele1_scPhi);
  ntu_MC->SetBranchAddress("ele1_EOverP",&ele1_EoP);
  ntu_MC->SetBranchAddress("ele1_tkP",   &ele1_tkP);
  ntu_MC->SetBranchAddress("ele1_e3x3",  &ele1_scE3x3);
  ntu_MC->SetBranchAddress("ele1_scERaw",&ele1_scERaw);
  ntu_MC->SetBranchAddress("ele1_scE",   &ele1_scE);
  ntu_MC->SetBranchAddress("ele1_scE_regression", &ele1_scE_regression);
  ntu_MC->SetBranchAddress("ele1_eSeedBC", &ele1_eSeedBC);
  ntu_MC->SetBranchAddress("ele1_charge",    &ele1_charge);
  ntu_MC->SetBranchAddress("ele1_scLocalPhi",&ele1_scLocalPhi); 
  ntu_MC->SetBranchAddress("ele1_scLocalEta",&ele1_scLocalEta); 
  ntu_MC->SetBranchAddress("ele1_bcN",&ele1_bcN); 
  if (useZ_){
    ntu_MC->SetBranchAddress("ele2_scEta", &ele2_scEta);
    ntu_MC->SetBranchAddress("ele2_scPhi", &ele2_scPhi);
    ntu_MC->SetBranchAddress("ele2_EOverP",&ele2_EoP);
    ntu_MC->SetBranchAddress("ele2_tkP",   &ele2_tkP);
    ntu_MC->SetBranchAddress("ele2_e3x3",  &ele2_scE3x3);
    ntu_MC->SetBranchAddress("ele2_scE",   &ele2_scE);
    ntu_MC->SetBranchAddress("ele2_scERaw",&ele2_scERaw);
    ntu_MC->SetBranchAddress("ele2_scE_regression", &ele2_scE_regression);
    ntu_MC->SetBranchAddress("ele2_eSeedBC", &ele2_eSeedBC);
    ntu_MC->SetBranchAddress("ele2_charge",    &ele2_charge);
    ntu_MC->SetBranchAddress("ele2_scLocalPhi",&ele2_scLocalPhi); 
    ntu_MC->SetBranchAddress("ele2_scLocalEta",&ele2_scLocalEta); 
    ntu_MC->SetBranchAddress("ele2_bcN",&ele2_bcN); 
  }

  //---- Set branch addresses for Data
  ActivateBranches(ntu_Data,branchesToRead,false);
  //ntu_Data->SetBranchAddress("PV_z", &PV_z);
  //ntu_Data->SetBranchAddress("PV_n", &PV_n);
  ntu_Data->SetBranchAddress("isZ", &isZ);
  ntu_Data->SetBranchAddress("ele1_scEta", &ele1_scEta);
  ntu_Data->SetBranchAddress("ele1_scPhi", &ele1_scPhi);
  ntu_Data->SetBranchAddress("ele1_EOverP",&ele1_EoP);
  ntu_Data->SetBranchAddress("ele1_tkP",   &ele1_tkP);
  ntu_Data->SetBranchAddress("ele1_eSeedBC",  &ele1_eSeedBC);
  ntu_Data->SetBranchAddress("ele1_e3x3",  &ele1_scE3x3);
  ntu_Data->SetBranchAddress("ele1_scERaw",&ele1_scERaw);
  ntu_Data->SetBranchAddress("ele1_scE",   &ele1_scE);
  ntu_Data->SetBranchAddress("ele1_scE_regression", &ele1_scE_regression);
  ntu_Data->SetBranchAddress("ele1_charge",    &ele1_charge);
  ntu_Data->SetBranchAddress("ele1_scLocalPhi",&ele1_scLocalPhi); 
  ntu_Data->SetBranchAddress("ele1_scLocalEta",&ele1_scLocalEta); 
  ntu_Data->SetBranchAddress("ele1_bcN",&ele1_bcN); 
  if (useZ_){
    ntu_Data->SetBranchAddress("ele2_scEta", &ele2_scEta);
    ntu_Data->SetBranchAddress("ele2_scPhi", &ele2_scPhi);
    ntu_Data->SetBranchAddress("ele2_EOverP",&ele2_EoP);
    ntu_Data->SetBranchAddress("ele2_tkP",   &ele2_tkP);
    ntu_Data->SetBranchAddress("ele2_eSeedBC",  &ele2_eSeedBC);
    ntu_Data->SetBranchAddress("ele2_e3x3",  &ele2_scE3x3);
    ntu_Data->SetBranchAddress("ele2_scERaw",&ele2_scERaw);
    ntu_Data->SetBranchAddress("ele2_scE",   &ele2_scE);
    ntu_Data->SetBranchAddress("ele2_scE_regression", &ele2_scE_regression);
    ntu_Data->SetBranchAddress("ele2_charge",    &ele2_charge);
    ntu_Data->SetBranchAddress("ele2_scLocalPhi",&ele2_scLocalPhi); 
    ntu_Data->SetBranchAddress("ele2_scLocalEta",&ele2_scLocalEta); 
    ntu_Data->SetBranchAddress("ele2_bcN",&ele2_bcN); 
  }
  const unsigned int nBins = 20;
  const int Ntempl = 4;
  std::cout << "nBins = " << nBins << std::endl;
  
  // histogram definition
  TH1F* h_EoP_MC[nBins][Ntempl] ;   
  TH1F* h_EoC_MC[nBins][Ntempl] ;
  TH1F* h_EoP_Data[nBins][Ntempl];
  TH1F* h_EoC_Data[nBins][Ntempl] ;
  
  for(int mod=0; mod<Ntempl; mod++){
    for(unsigned int i = 0; i < nBins; ++i)
      {
	char histoName[80];
	sprintf(histoName, "EoP_MC_%d_mod%d", i,mod+1);
	h_EoP_MC[i][mod] = new TH1F(histoName, histoName, 1200, 0., 3.);
	h_EoP_MC[i][mod] -> SetFillColor(4);
	h_EoP_MC[i][mod] -> SetFillStyle(3004);
	sprintf(histoName, "EoC_MC_%d_mod%d", i,mod+1);
	h_EoC_MC[i][mod] = new TH1F(histoName, histoName, 1200, 0., 3.);
	h_EoC_MC[i][mod] -> SetFillColor(3);
	h_EoC_MC[i][mod] -> SetFillStyle(3004);
	sprintf(histoName, "EoP_Data_%d_mod%d", i,mod+1);
	h_EoP_Data[i][mod] = new TH1F(histoName, histoName, 1200, 0., 3.);
	h_EoP_Data[i][mod] -> SetFillColor(4);
	h_EoP_Data[i][mod] -> SetFillStyle(3004);
	sprintf(histoName, "EoC_Data_%d_mod%d", i,mod+1);
	h_EoC_Data[i][mod] = new TH1F(histoName, histoName, 1200, 0., 3.);
	h_EoC_Data[i][mod] -> SetFillColor(3);
	h_EoC_Data[i][mod] -> SetFillStyle(3004);
      }
  }

  TH1F* h_spreadEoP_MC[Ntempl] ;   
  TH1F* h_spreadEoC_MC[Ntempl] ;
  TH1F* h_spreadEoP_Data[Ntempl] ;   
  TH1F* h_spreadEoC_Data[Ntempl] ;
  TH1F* h_spreadRatio[Ntempl] ;
  TH1F* h_spreadRatioCorr[Ntempl] ;

  for(int mod=0; mod<Ntempl; mod++){
    char histoName[80];
    sprintf(histoName, "spreadEoP_MC_mod%d",mod+1);
    h_spreadEoP_MC[mod] = new TH1F(histoName, histoName, 400, 0.9, 1.1);
    h_spreadEoP_MC[mod] -> SetFillColor(4);
    h_spreadEoP_MC[mod] -> SetFillStyle(3004);

    sprintf(histoName, "spreadEoC_MC_mod%d",mod+1);
    h_spreadEoC_MC[mod] = new TH1F(histoName, histoName, 400, 0.9, 1.1);
    h_spreadEoC_MC[mod] -> SetFillColor(3);
    h_spreadEoC_MC[mod] -> SetFillStyle(3004);

    sprintf(histoName, "spreadEoP_Data_mod%d",mod+1);
    h_spreadEoP_Data[mod] = new TH1F(histoName, histoName, 400, 0.9, 1.1);
    h_spreadEoP_Data[mod] -> SetFillColor(4);
    h_spreadEoP_Data[mod] -> SetFillStyle(3004);

    sprintf(histoName, "spreadEoC_Data_mod%d", mod+1);
    h_spreadEoC_Data[mod] = new TH1F(histoName, histoName, 400, 0.9, 1.1);
    h_spreadEoC_Data[mod] -> SetFillColor(3);
    h_spreadEoC_Data[mod] -> SetFillStyle(3004);

    sprintf(histoName, "spreadRatio_mod%d", mod+1);
    h_spreadRatio[mod] = new TH1F(histoName, histoName, 400, 0.9, 1.1);
    h_spreadRatio[mod] -> SetFillColor(3);
    h_spreadRatio[mod] -> SetFillStyle(3004);

    sprintf(histoName, "spreadRatioCorr_mod%d", mod+1);
    h_spreadRatioCorr[mod] = new TH1F(histoName, histoName, 400, 0.9, 1.1);
    h_spreadRatioCorr[mod] -> SetFillColor(3);
    h_spreadRatioCorr[mod] -> SetFillStyle(3004);

  }


  //---- book templates
  // [0] --> uncorrected
  // [1] --> corrected
  TH1F* h_template_MC[Ntempl][2];
  TH1F* h_template_Data[Ntempl][2];
  for(unsigned int i = 0; i < Ntempl; ++i){
    char histoName[100];
    sprintf(histoName, "template_MC_%d", i);
    h_template_MC[i][0] = new TH1F(histoName, "", 1200, 0., 3.);   
    sprintf(histoName, "template_DATA_%d", i);
    h_template_Data[i][0] = new TH1F(histoName, "", 1200, 0., 3.);   
    sprintf(histoName, "template_MCcorr_%d", i);
    h_template_MC[i][1] = new TH1F(histoName, "", 1200, 0., 3.);   
    sprintf(histoName, "template_DATAcorr_%d", i);
    h_template_Data[i][1] = new TH1F(histoName, "", 1200, 0., 3.);   
  }


  TH1F *hnbc = new TH1F ("hnbc","hnbc",10,0,10);

  //******************************************************************************************
  //*************************************** MC  ********************************************** 
  std::cout << "Loop on MC events ... " << std::endl; 
 
  float ww = 1 ;
  //---- loop on MC, make refernce and fit dist
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    
    if( entry%500000 == 0 ) std::cout << "reading saved entry " << entry << std::endl;
    //if (entry>200000) break;
    ntu_MC->GetEntry(entry);

    // -- PU weights
    if (usePUweights_) ww = puReweighting->GetWeight(npu);

    //-- eta or R9 cuts
    ele1_R9 = ele1_scE3x3/ele1_scERaw;
 
    // --- first electron
    
    //-- phi cracks
    float phi = (ele1_scPhi+PI)/XTALWIDTH;;
    float modphi = (int)phi%20;

    // -- eta gaps
    float fetaCry = fabs (ele1_scEta) / XTALWIDTH;
    
    //-- local eta 
    float locEta = ele1_scLocalEta+0.5; 
    if (useLocalPhi_) locEta = ele1_scLocalPhi+0.5;
    
    int mod, bin;
   
    if ( fabs(ele1_scEta) < etaMax_ && (ele1_R9 > r9min_ && ele1_R9 < r9max_ ) &&  (ele1_bcN >= bcNmin_ && ele1_bcN <= bcNmax_) &&
	 fabs(modphi-10)>2 && IsEtaGap(fetaCry) == 0 && fabs(locEta-0.5) < 0.5) {
      
      mod = templIndex(fetaCry);
      bin = nBins * (locEta); 
      //FillHistograms(ele1_scE, ele1_scE_regression, ele1_tkP, h_template_MC[mod][0], h_template_MC[mod][1], h_EoP_MC[bin][mod], h_EoC_MC[bin][mod], ww, mod, bin, nBins);
      //FillHistograms(ele1_EoP*ele1_tkP, ele1_scE_regression, ele1_tkP, h_template_MC[mod][0], h_template_MC[mod][1], h_EoP_MC[bin][mod], h_EoC_MC[bin][mod], ww, mod, bin, nBins);
      FillHistograms(ele1_EoP*ele1_tkP, ele1_scE, ele1_tkP, h_template_MC[mod][0], h_template_MC[mod][1], h_EoP_MC[bin][mod], h_EoC_MC[bin][mod], ww, mod, bin, nBins);
    }

     
    // --- second  electron
    /*  
    if (isZ) {
    
      ele2_R9 = ele2_scE3x3/ele2_scERaw;

      
      //-- phi cracks
      phi = (ele2_scPhi+3.1415926536)/XTALWIDTH;;
      modphi = (int)phi%20;
      
      // -- eta gaps
      fetaCry = fabs (ele2_scEta) / XTALWIDTH;
      
      //-- local eta 
      locEta = ele2_scLocalEta+0.5; 
      if (useLocalPhi) locEta = ele2_scLocalPhi+0.5;

      if ( fabs(ele2_scEta) < etaMax && (ele2_R9 > r9min && ele2_R9 < r9max ) && (ele2_bcN >= bcNmin && ele2_bcN <= bcNmax) &&
	   fabs(modphi-10)>2 && IsEtaGap(fetaCry) == 0 && fabs(locEta-0.5) < 0.5) {
          
	//-- fill templates for each mod
	mod = templIndex(fetaCry);
	correction = ele2_scE_regression/ele2_scE;
	h_template_MC[mod][0]-> Fill(ele2_EoP,ww);
	h_template_MC[mod][1]-> Fill(ele2_EoP*correction,ww);
	
	//-- fill MC histos in eta bins
	bin = nBins * (locEta); 
	if (bin>nBins-1 || bin < 0 ) {
	  std::cout << "Error in bins: " << bin << " " << ele2_scLocalEta << std::endl;
	  continue;
	}
	
	h_EoP_MC[bin][mod] -> Fill(ele2_EoP,ww);
	h_EoC_MC[bin][mod] -> Fill(ele2_EoP*correction,ww);
      }
    } // end second ele
    */

  }

  std::cout << "stat: "<< h_EoP_MC[0][10]->GetEntries()<< " " <<  h_EoC_MC[0][10]->GetEntries() << std::endl; 

  delete   ntu_MC;

  //******************************************************************************************
  //*************************************** DATA ********************************************** 
  std::cout << "Loop on Data events ..." << std::endl; 
  //--- loop on data
  for(int entry = 0; entry < ntu_Data->GetEntries(); ++entry) {
    if( entry%500000 == 0 ) std::cout << "reading saved entry " << entry << std::endl;
    //if (entry>200000) break;
    //if (entry%100!=0) continue;
    ntu_Data->GetEntry(entry);


    //if (fabs(PV_z)>3) continue;
    //if (PV_n<20) continue;

    ww= 1.;

    //-- eta or R9 cuts
    ele1_R9 = ele1_scE3x3/ele1_scERaw;


    // --- first electron
    
    //-- phi cracks
    float phi = (ele1_scPhi+PI)/XTALWIDTH;;
    float modphi = (int)phi%20;

    // -- eta gaps
    float fetaCry = fabs (ele1_scEta) / XTALWIDTH;

    //-- local eta 
    float locEta = ele1_scLocalEta+0.5; 
    if (useLocalPhi_) locEta = ele1_scLocalPhi+0.5;

    int mod, bin;

    if ( fabs(ele1_scEta) < etaMax_ && (ele1_R9 > r9min_ && ele1_R9 < r9max_ ) && (ele1_bcN >= bcNmin_ && ele1_bcN <= bcNmax_) &&
	 fabs(modphi-10)>2 && IsEtaGap(fetaCry) == 0 && fabs(locEta-0.5) < 0.5) {
      
      hnbc ->Fill(ele1_bcN );

      //-- fill templates for each mod
      mod = templIndex(fetaCry);
      //if ( applyLocalCorrections_) {
      //	//float newE = ele1_scE - ele1_eSeedBC + ele1_eSeedBC/( feta[mod]->Eval(ele1_scLocalEta)* fphi[mod]->Eval(ele1_scLocalPhi));
      //	//float newE = ele1_scE - ele1_eSeedBC + ele1_eSeedBC/( feta[mod]->Eval(ele1_scLocalEta));
      //	float newE = ele1_scE - ele1_eSeedBC + ele1_eSeedBC/( fphi[mod]->Eval(ele1_scLocalPhi));
      //	ele1_EoP = ele1_EoP/ele1_scE * newE;
      //      }

      bin = nBins * (locEta); 
      //FillHistograms(ele1_scE, ele1_scE_regression, ele1_tkP, h_template_Data[mod][0], h_template_Data[mod][1], h_EoP_Data[bin][mod], h_EoC_Data[bin][mod], ww, mod, bin, nBins);
      // in 53X rereco scE is already corrected for local containment, EoP is not --> use EoP to look at dependence on local coordinate
      //FillHistograms(ele1_EoP*ele1_tkP, ele1_scE_regression, ele1_tkP, h_template_Data[mod][0], h_template_Data[mod][1], h_EoP_Data[bin][mod], h_EoC_Data[bin][mod], ww, mod, bin, nBins);
      FillHistograms(ele1_EoP*ele1_tkP, ele1_scE, ele1_tkP, h_template_Data[mod][0], h_template_Data[mod][1], h_EoP_Data[bin][mod], h_EoC_Data[bin][mod], ww, mod, bin, nBins);

    }

     
    // --- second  electron
    /*
    if (isZ) {
      ele2_R9 = ele2_scE3x3/ele2_scERaw;
      
      //-- phi cracks
      phi = (ele2_scPhi+3.1415926536)/XTALWIDTH;;
      modphi = (int)phi%20;

      // -- eta gaps
      fetaCry = fabs (ele2_scEta) / XTALWIDTH;

      //-- local eta 
      locEta = ele2_scLocalEta+0.5; 
      if (useLocalPhi) locEta = ele1_scLocalPhi+0.5;

      if ( fabs(ele2_scEta) < etaMax && (ele2_R9 > r9min && ele2_R9 < r9max ) && (ele2_bcN >= bcNmin && ele2_bcN <= bcNmax) &&
	   fabs(modphi-10)>2 && IsEtaGap(fetaCry) == 0 && fabs(locEta-0.5) < 0.5) {

	//-- fill templates for each mod
	mod = templIndex(fetaCry);
	correction = ele2_scE_regression/ele2_scE;
	h_template_Data[mod][0]-> Fill(ele2_EoP,ww);
	h_template_Data[mod][1]-> Fill(ele2_EoP*correction,ww);
	
	//-- fill MC histos in eta bins
	bin = nBins * (locEta); 
	if (bin>nBins-1 || bin < 0 ) {
	  std::cout << "Error in bins: " << bin << " " << ele2_scLocalEta << std::endl;
	  continue;
	}
	
	h_EoP_Data[bin][mod] -> Fill(ele2_EoP,ww);
	h_EoC_Data[bin][mod] -> Fill(ele2_EoP*correction,ww);
      }
    } // end second ele

    */
  }
  

  delete   ntu_Data;
  
  ///////////////****************** Fit the histograms and fill the graphs *************** ////////////////////////
  int rebin = 4;

  TGraphErrors* g_EoP_MC[Ntempl];
  TGraphErrors* g_EoC_MC[Ntempl];
  TGraphErrors* g_ratio_MC[Ntempl];
  
  TGraphErrors* g_EoP_Data[Ntempl];
  TGraphErrors* g_EoC_Data[Ntempl];
  TGraphErrors* g_ratio_Data[Ntempl];

  TGraphErrors* g_ratio_uncorr[Ntempl];
  TGraphErrors* g_ratio_corr[Ntempl];
 
  for (int mod=0; mod<4; mod++){
    char histoName[100];
    
    sprintf(histoName, "gEoP_MC_mod%d", mod+1);
    g_EoP_MC[mod]   = new TGraphErrors(); 
    g_EoP_MC[mod]->SetName(histoName);
    g_EoP_MC[mod]->SetLineColor(kRed);
    
    sprintf(histoName, "gEoC_MC_mod%d", mod+1);
    g_EoC_MC[mod]   = new TGraphErrors(); 
    g_EoC_MC[mod]->SetName(histoName);
    g_EoC_MC[mod]->SetLineColor(kRed+2);
    
    sprintf(histoName, "gRatio_MC_mod%d", mod+1);
    g_ratio_MC[mod]   = new TGraphErrors(); 
    g_ratio_MC[mod]->SetName(histoName);
    
    sprintf(histoName, "gEoP_Data_mod%d", mod+1);
    g_EoP_Data[mod]   = new TGraphErrors(); 
    g_EoP_Data[mod]->SetName(histoName);
    g_EoP_Data[mod]->SetLineColor(kGreen);
    
    sprintf(histoName, "gEoC_Data_mod%d", mod+1);
    g_EoC_Data[mod]   = new TGraphErrors(); 
    g_EoC_Data[mod]->SetName(histoName);
    g_EoC_Data[mod]->SetLineColor(kGreen+2);

    sprintf(histoName, "gRatio_Data_mod%d", mod+1);
    g_ratio_Data[mod]   = new TGraphErrors(); 
    g_ratio_Data[mod]->SetName(histoName);

    sprintf(histoName, "gRatio_uncorr_mod%d", mod+1);
    g_ratio_uncorr[mod]   = new TGraphErrors(); 
    g_ratio_uncorr[mod]->SetName(histoName);

    sprintf(histoName, "gRatio_corr_mod%d", mod+1);
    g_ratio_corr[mod]   = new TGraphErrors(); 
    g_ratio_corr[mod]->SetName(histoName);
 
  }
  
  

  //************************************* FITTING ***************************************************//
  
  float eopmin = 0.7;
  float eopmax = 1.5;
  
  TF1 *fitfun = new TF1("fitfun", crystalBallLowHigh, eopmin, eopmax, 7);
  fitfun -> SetNpx(10000);
  fitfun -> SetLineWidth(2);
  fitfun -> SetLineColor(kRed);
  //  fitfun -> FixParameter(7,1.);
  fitfun -> SetParLimits(3,0.,10.);
  fitfun -> SetParLimits(4,0.,10.);
  fitfun -> SetParLimits(5,0.,10.);
  
  fitfun -> SetParName(0,"N");
  fitfun -> SetParName(1,"#mu");
  fitfun -> SetParName(2,"#sigma");
  fitfun -> SetParName(3,"#alpha_{high}");
  fitfun -> SetParName(4,"n_{high}");
  fitfun -> SetParName(5,"#alpha_{low}");
  fitfun -> SetParName(6,"n_{low}");

  float initialPars[7]={100.,1.,0.05,0.8,2.,1.5,3.5};

  for(int mod=0;mod<4;mod++){
    
    for(unsigned int i = 0; i < nBins; ++i)
      {
	std::cout << "***** Fitting :  mod: " <<mod <<"  bin:  "<< i << std::endl;
	
	h_EoP_MC[i][mod] -> Rebin(rebin);    
	h_EoC_MC[i][mod] -> Rebin(rebin);    
	h_EoP_Data[i][mod] -> Rebin(rebin);    
	h_EoC_Data[i][mod] -> Rebin(rebin);    
	
       	float xval = (i+0.5)*1/(float)nBins - 0.5;

	//************************ MC ****************************************************************
	
	//--- uncorrected MC
	float scaleMC, escaleMC;
	DoFit(h_EoP_MC[i][mod], fitfun, scaleMC, escaleMC);
      	g_EoP_MC[mod] -> SetPoint(i, xval , scaleMC );
	g_EoP_MC[mod] -> SetPointError(i, 0., escaleMC );
	h_spreadEoP_MC[mod] -> Fill(scaleMC);

	//--- corrected MC   
	float scaleMCcorr, escaleMCcorr;
	DoFit(h_EoC_MC[i][mod], fitfun, scaleMCcorr, escaleMCcorr);
      	g_EoC_MC[mod] -> SetPoint(i, xval , scaleMCcorr );
	g_EoC_MC[mod] -> SetPointError(i, 0., escaleMCcorr );
	h_spreadEoC_MC[mod] -> Fill(scaleMCcorr);

	//************************ DATA **********************************************
	//--- uncorrected data
	float scaleDA, escaleDA;
	DoFit(h_EoP_Data[i][mod], fitfun, scaleDA, escaleDA);
      	g_EoP_Data[mod] -> SetPoint(i, xval , scaleDA );
	g_EoP_Data[mod] -> SetPointError(i, 0., escaleDA );
	h_spreadEoP_Data[mod] -> Fill(scaleDA);

	//--- corrected Data 
	float scaleDAcorr, escaleDAcorr;
	DoFit(h_EoC_Data[i][mod], fitfun, scaleDAcorr, escaleDAcorr);
      	g_EoC_Data[mod] -> SetPoint(i, xval , scaleDAcorr );
	g_EoC_Data[mod] -> SetPointError(i, 0., escaleDAcorr );
	h_spreadEoC_Data[mod] -> Fill(scaleDAcorr);

	//--- ratio finalization MC corr/uncorr
	float ratioMC = scaleMCcorr/scaleMC;
	float eratioMC = ratioMC*sqrt(pow(escaleMC/scaleMC,2) + pow(escaleMCcorr/scaleMCcorr,2)); 
	g_ratio_MC[mod] -> SetPoint(i, xval, ratioMC);
	g_ratio_MC[mod] -> SetPointError(i, 0., eratioMC);

	//--- ratio finalization DATA corr/uncorr
	float ratioDA  = scaleDAcorr/scaleDA;
	float eratioDA = ratioDA*sqrt(pow(escaleDA/scaleDA,2) + pow(escaleDAcorr/scaleDAcorr,2)); 
	g_ratio_Data[mod] -> SetPoint(i, xval , ratioDA);
	g_ratio_Data[mod] -> SetPointError(i, 0., eratioDA);


	//--- ratio finalization data/MC uncorrected
	float ratioU = scaleDA/scaleMC;
	float eratioU = ratioU*sqrt(pow(escaleMC/scaleMC,2) + pow(escaleDA/scaleDA,2)); 
	g_ratio_uncorr[mod] -> SetPoint(i, xval , ratioU);
	g_ratio_uncorr[mod] -> SetPointError(i, 0., eratioU);
	h_spreadRatio[mod] -> Fill(ratioU);
    

	//--- ratio finalization data/MC corrected
	float ratioC = scaleDAcorr/scaleMCcorr;
	float eratioC = ratioC*sqrt(pow(escaleMCcorr/scaleMCcorr,2) + pow(escaleDAcorr/scaleDAcorr,2)); 
	g_ratio_corr[mod] -> SetPoint(i, xval , ratioC);
	g_ratio_corr[mod] -> SetPointError(i, 0., eratioC);
    	h_spreadRatioCorr[mod] -> Fill(ratioC);
      }
  }

  
  TFile fout(outfilename_.c_str(),"recreate");
  hnbc->Write();
  for(int mod=0;mod<4;mod++) {
    g_EoP_MC[mod]->Write();
    g_EoC_MC[mod]->Write();
    g_EoP_Data[mod]->Write();
    g_EoC_Data[mod]->Write();
    g_ratio_MC[mod]->Write();
    g_ratio_Data[mod]->Write();
    g_ratio_uncorr[mod]->Write();
    g_ratio_corr[mod]->Write();

    h_spreadEoP_MC[mod]->Write();
    h_spreadEoC_MC[mod]->Write();
    h_spreadEoP_Data[mod]->Write();
    h_spreadEoC_Data[mod]->Write();

    h_spreadRatio[mod]->Write();
    h_spreadRatioCorr[mod]->Write();
 
    h_template_MC[mod][0]-> Write();
    h_template_MC[mod][1]-> Write();
    h_template_Data[mod][0]-> Write();
    h_template_Data[mod][1]-> Write();
    for (int ibin = 0; ibin < nBins ; ibin ++){
      h_EoP_MC[ibin][mod] ->Write();
      h_EoP_Data[ibin][mod] ->Write();
      h_EoC_MC[ibin][mod] ->Write();
      h_EoC_Data[ibin][mod] ->Write();
    }
  }
  fout.Close();
}
