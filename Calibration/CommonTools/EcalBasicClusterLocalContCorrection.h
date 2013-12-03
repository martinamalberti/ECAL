#ifndef RecoEcal_EgammaCoreTools_EcalBasicClusterLocalContCorrection_h
#define RecoEcal_EgammaCoreTools_EcalBasicClusterLocalContCorrection_h

/** \class EcalBasicClusterLocalContCorrection
  *  Function to correct basic cluster energy for local containment
  *
  *  $Id: EcalBasicClusterLocalContCorrection.h
  *  $Date:
  *  $Revision:
  *  \author Martina Malberti, February 2012
  */

#include "RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrectionBaseClass.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"


class EcalBasicClusterLocalContCorrection : public EcalClusterLocalContCorrectionBaseClass {
 
 public:
  EcalBasicClusterLocalContCorrection( const edm::ParameterSet &) {};
  
  // compute the correction
  float getCorrectionVsLocalEta( const CaloGeometry *caloGeometry ,reco::CaloCluster & ) ;
  float getCorrectionVsLocalPhi( const CaloGeometry *caloGeometry ,reco::CaloCluster & ) ;

 private:
  int   getEcalModule( DetId id );
  std::pair<double,double> getLocalPosition(const CaloGeometry *caloGeometry, reco::CaloCluster & caloCluster) ;
};

#endif
