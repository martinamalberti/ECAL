// Local containment correction vs local eta 
// --- correction function derived from E/p data on Zee 2012
// --- tested both pol2
  
#include "TF1.h"

//TF1 *f = new TF1("f","[0]+([1]*(x)+[2]*pow(x,2))",-1,1.)
TF1 *f = new TF1(f,"[0]+([1]*(x-0.0)+[2]*pow(x-0.0,2)+[3]*pow(x-0.0,3))",-1,1.);

double LocalContainmentCorrectionEta (float localEta, int imod, bool data){
  


  // pol 2
//   // mc corr
//   if (imod==0) f->SetParameters(1.00197,0.000344159 , -0.0197935);
//   if (imod==1) f->SetParameters(1.001,-0.000985039 , -0.00654623);
//   if (imod==2) f->SetParameters(1.0008,-0.000436675 , -0.00891791);
//   if (imod==3) f->SetParameters(1.000, 0, 0);

//   // data corr
//   if (data){
//     if (imod==0) f->SetParameters(1.00353,0.00287393 , -0.0411409);
//     if (imod==1) f->SetParameters(1.00346,0.00286825 , -0.0418321);
//     if (imod==2) f->SetParameters(1.00242,0.00724232 , -0.0371932);
//     if (imod==3) f->SetParameters(1.0,0.0, 0.0);
//   }


  // pol3
  // mc corr
  if (imod==0)f->SetParameters(1.00184,-0.000679226 , -0.0205473 , 0.00193425);
  if (imod==1)f->SetParameters(1.00111,0.00155065 , -0.0116107 , -0.00971859);
  if (imod==2)f->SetParameters(1.00043,0.000776506 , -0.0105513 , 0.00204469);
  if (imod==3)f->SetParameters(1.0,0. , 0. , 0);

  // data corr
  if (data){
    if (imod==0) f->SetParameters(1.00353,0.00931611 , -0.0409904 , -0.0413427);
    if (imod==2) f->SetParameters(1.00344,0.00783796 , -0.0427353 , -0.0344312);
    if (imod==3) f->SetParameters(1.00376,0.00932689 , -0.0509447 , -0.0182674);
    if (imod==4) f->SetParameters(1.0, 0. , 0.0 , 0.0);
  }


  double corr = f-> Eval(localEta);
  return(1./corr);

}

