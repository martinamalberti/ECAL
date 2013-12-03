#include "templateUtils.h"

/*** double crystall ball ***/
double crystalBallLowHigh(double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  //[5] = alpha2
  //[6] = n2
  //[7] = scale
  
  //double xx = x[0] * par[7];
  double xx = x[0];
  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];
  double n = par[4];
  double alpha2 = par[5];
  double n2 = par[6];
  
  if( (xx-mean)/sigma > fabs(alpha) )  
  {
    double A = pow(n/fabs(alpha), n) * exp(-0.5 * alpha*alpha);
    double B = n/fabs(alpha) - fabs(alpha);
    
    //    return par[0] * par[7] * A * pow(B + (xx-mean)/sigma, -1.*n);
    return par[0] * A * pow(B + (xx-mean)/sigma, -1.*n);
  }
  
  else if( (xx-mean)/sigma < -1.*fabs(alpha2) )
  {
    double A = pow(n2/fabs(alpha2), n2) * exp(-0.5 * alpha2*alpha2);
    double B = n2/fabs(alpha2) - fabs(alpha2);
    
    //    return par[0] * par[7] * A * pow(B - (xx-mean)/sigma, -1.*n2);
    return par[0] * A * pow(B - (xx-mean)/sigma, -1.*n2);
  }
  
  else
  {
    //return par[0] * par[7] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigma*sigma) );
    return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigma*sigma) );
  } 
  
}

