// Local containment correction vs local phi
// --- correction function derived from E/p data on Zee 2012
// --- tested both pol2
  
#include "TF1.h"

//TF1 *f = new TF1("f","[0]+([1]*(x)+[2]*pow(x,2))",-1,1.)
TF1 *f = new TF1(f,"[0]+([1]*(x-0.0)+[2]*pow(x-0.0,2)+[3]*pow(x-0.0,3))",-1,1.);

double LocalContainmentCorrectionPhi (float localPhi, int imod, bool data){
  

  //MC
  if (imod==0) f->SetParameters(1.00239,0.00058966 , -0.0270848 , -0.00253261);
  if (imod==1) f->SetParameters(1.00351,-0.00246757 , -0.0441807 , 0.0093244);
  if (imod==2) f->SetParameters(1.00283,-0.00459891 , -0.0335485 , 0.0269317);
  if (imod==3) f->SetParameters(1,0. , 0.0 , 0.);
  
  //DATA
  if (data){
    if (imod==0) f->SetParameters(1.00302,-0.00260969 , -0.0336908 , 0.0182451);
    if (imod==1) f->SetParameters(1.00339,-0.0010086 , -0.0384517 , 0.0100485);
    if (imod==2) f->SetParameters(1.00349,-0.000234313 , -0.042849 , 0.00813363);
    if (imod==3) f->SetParameters(1,0. , 0.0 , 0.);
  }


  double corr = f-> Eval(localPhi);
  return(1./corr);

}

