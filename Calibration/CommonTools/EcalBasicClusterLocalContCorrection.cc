#include "RecoEcal/EgammaCoreTools/plugins/EcalBasicClusterLocalContCorrection.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"

#include "TVector2.h"
#include "TMath.h"
#include "TF1.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <iostream>

using namespace std;
using namespace edm;


//****** Containment correction vs local eta
float EcalBasicClusterLocalContCorrection::getCorrectionVsLocalEta( const CaloGeometry *caloGeometry , reco::CaloCluster & basicCluster ) 
{

  //--- local cluster coordinates 
  std::pair<double,double> localPosition = getLocalPosition(caloGeometry, basicCluster);
  float localEta = localPosition.first;
  
  //--- ecal module
  int imod    = getEcalModule(basicCluster.seed());

  //--- data driven correction from E/p (parametrized with pol2)
  TF1 *f = new TF1("f","[0]+([1]*(x)-[2]*pow(x,2))",-1,1.);
  if (imod==0) f->SetParameters(1.00603,0.00300789 , 0.0667232);
  if (imod==1) f->SetParameters(1.00655,0.00386189 , 0.073931);
  if (imod==2) f->SetParameters(1.00634,0.00631341 , 0.0764134);
  if (imod==3) f->SetParameters(1.00957,0.0113306 , 0.123808);

  double corr = f-> Eval(localEta);
  return(1./corr);
 
}

//****** Containment correction vs local phi
float EcalBasicClusterLocalContCorrection::getCorrectionVsLocalPhi( const CaloGeometry *caloGeometry , reco::CaloCluster & basicCluster ) 
{

  //--- local cluster coordinates 
  std::pair<double,double> localPosition = getLocalPosition(caloGeometry, basicCluster);
  float localPhi = localPosition.second;
  
  //--- global cluster coordinates 
  int imod    = getEcalModule(basicCluster.seed());

  //--- data driven correction from E/p (parametrized with pol2)
  TF1 *f = new TF1("f","[0]+([1]*(x)-[2]*pow(x,2))",-1,1.);
  
  if (imod==0) f->SetParameters(1.00403,-0.0012733 , 0.042925);
  if (imod==1) f->SetParameters(1.00394,-0.00137567 , 0.0416698);
  if (imod==2) f->SetParameters(1.00298,-0.00111589 , 0.0320377);
  if (imod==3)f->SetParameters(1.00269,-0.00153347 , 0.0296769);

  double corr = f-> Eval(localPhi);
  return(1./corr);
 
}

//****** Cluster position in ECAL
int EcalBasicClusterLocalContCorrection::getEcalModule( DetId id )
{
  int mod = 0;

  //  int ieta = (EBDetId(id.first)).ieta();
  int ieta = (EBDetId(id)).ieta();

  if (fabs(ieta) <=25 ) mod = 0;
  if (fabs(ieta) > 25 && fabs(ieta) <=45 ) mod = 1;
  if (fabs(ieta) > 45 && fabs(ieta) <=65 ) mod = 2;
  if (fabs(ieta) > 65 && fabs(ieta) <=85 ) mod = 3;

  return (mod);
 
}




//****** Cluster local coordinates
std::pair<double,double> EcalBasicClusterLocalContCorrection::getLocalPosition(const CaloGeometry *caloGeometry, reco::CaloCluster & caloCluster)

{
  //--------------if barrel calculate local position wrt xtal center -------------------
  const CaloSubdetectorGeometry* geom = caloGeometry->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);//EcalBarrel = 1
  
  const math::XYZPoint position_ = caloCluster.position(); 
  double Theta = -position_.theta()+0.5*TMath::Pi();
  double Eta = position_.eta();
  double Phi = TVector2::Phi_mpi_pi(position_.phi());
  
  //Calculate expected depth of the maximum shower from energy (like in PositionCalc::Calculate_Location()):
  // The parameters X0 and T0 are hardcoded here because these values were used to calculate the corrections:
  const float X0 = 0.89; const float T0 = 7.4;
  double depth = X0 * (T0 + log(caloCluster.energy()));
  
  
  //search which crystal is closest to the cluster position and call it crystalseed:
  //std::vector<DetId> crystals_vector = *scRef.getHitsByDetId();   //deprecated
  std::vector< std::pair<DetId, float> > crystals_vector = caloCluster.hitsAndFractions();
  float dphimin=999.;
  float detamin=999.;
  int ietaclosest = 0;
  int iphiclosest = 0;
  for (unsigned int icry=0; icry!=crystals_vector.size(); ++icry) 
    {    
      EBDetId crystal(crystals_vector[icry].first);
      const CaloCellGeometry* cell=geom->getGeometry(crystal);
      GlobalPoint center_pos = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(depth);
      double EtaCentr = center_pos.eta();
      double PhiCentr = TVector2::Phi_mpi_pi(center_pos.phi());
      if (TMath::Abs(EtaCentr-Eta) < detamin) {
	detamin = TMath::Abs(EtaCentr-Eta); 
	ietaclosest = crystal.ieta();
      }
      if (TMath::Abs(TVector2::Phi_mpi_pi(PhiCentr-Phi)) < dphimin) {
	dphimin = TMath::Abs(TVector2::Phi_mpi_pi(PhiCentr-Phi)); 
	iphiclosest = crystal.iphi();
      }
    }
  EBDetId crystalseed(ietaclosest, iphiclosest);
  
  // Get center cell position from shower depth
  const CaloCellGeometry* cell=geom->getGeometry(crystalseed);
  GlobalPoint center_pos = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(depth);
  
  //PHI
  double PhiCentr = TVector2::Phi_mpi_pi(center_pos.phi());
  double PhiWidth = (TMath::Pi()/180.);
  double PhiCry = (TVector2::Phi_mpi_pi(Phi-PhiCentr))/PhiWidth;
  if (PhiCry>0.5) PhiCry=0.5;
  if (PhiCry<-0.5) PhiCry=-0.5;
  //flip to take into account ECAL barrel symmetries:
  if (ietaclosest<0) PhiCry *= -1.;
  
  //ETA
  double ThetaCentr = -center_pos.theta()+0.5*TMath::Pi();
  double ThetaWidth = (TMath::Pi()/180.)*TMath::Cos(ThetaCentr);
  double EtaCry = (Theta-ThetaCentr)/ThetaWidth;    
  if (EtaCry>0.5) EtaCry=0.5;
  if (EtaCry<-0.5) EtaCry=-0.5;
  //flip to take into account ECAL barrel symmetries:
  if (ietaclosest<0) EtaCry *= -1.;
  
  //-------------- end calculate local position -------------
  
  std::pair<double,double> etaphi(EtaCry,PhiCry);
  return etaphi;
}

