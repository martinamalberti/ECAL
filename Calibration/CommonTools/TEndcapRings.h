#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TTree.h"
#include "TVirtualFitter.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TChain.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

class TEndcapRings {
 private:
  float iEndcapRing[100][100][2]; 
 public:
  TEndcapRings(); 
  ~TEndcapRings();
  Int_t GetEndcapRing(Int_t,Int_t,Int_t);
  Int_t GetEndcapIeta(Int_t,Int_t,Int_t);

  // ClassDef(TEndcapRings,1); //ring class
};

// default constructor, reading the map from file
TEndcapRings::TEndcapRings() {
  FILE *fRing;
  fRing = fopen("../CommonTools/eerings.dat","r");
  std::cout << "Inizializing endcap geometry from: eerings.dat" << std::endl;
  int ix,iy,iz,ir;
  while(fscanf(fRing,"(%d,%d,%d) %d \n",&ix,&iy,&iz,&ir) !=EOF ) {
    if (iz<0) iz=0; 
    iEndcapRing[ix][iy][iz] = ir;
  }
  return; 
}

TEndcapRings::~TEndcapRings() { return;}

Int_t TEndcapRings::GetEndcapRing(Int_t ix, Int_t iy, Int_t iz){
  return iEndcapRing[ix][iy][iz];
}

Int_t TEndcapRings::GetEndcapIeta(Int_t ix, Int_t iy, Int_t iz){
  Int_t iSide = iz; 
  if (iSide<0) iSide=0; 
  Int_t iEtaOffset = 86*iz; 
  Int_t iEta = iEtaOffset + iz*iEndcapRing[ix][iy][iSide];
  return iEta;
}
