//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Finite Deformation Hyper-Elastic classes
//# CLASS:             
//#
//# VERSION:           0.6_(1803398874989) (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              19AUg2003
//# UPDATE HISTORY:    Sept 2003
//#		       28May2004
//#
//===============================================================================

#ifndef SimoPisterWEnergy_CPP
#define SimoPisterWEnergy_CPP

#include <SimoPisterWEnergy.h>

//================================================================================
// Normal constructor
//================================================================================
SimoPisterWEnergy::SimoPisterWEnergy(double K_in) :K(K_in)
{

}

SimoPisterWEnergy::SimoPisterWEnergy( ) :K(0.0)
{

}

//================================================================================
// Normal destructor
//================================================================================
SimoPisterWEnergy::~SimoPisterWEnergy( )
{

}

//================================================================================
//create a clone of itself
//================================================================================
WEnergy * SimoPisterWEnergy::newObj( )
  {
    WEnergy  *new_WEnergy = new SimoPisterWEnergy(K);
    return new_WEnergy;
  }


//================================================================================
// w
//================================================================================
const double  SimoPisterWEnergy::wE( const double &J_in, const Vector &lambda_wave_in )
{
    double wEnergy = 0.25 * K * (J_in*J_in - 1.0 - 2.0*log(J_in));
    return wEnergy;
}

//================================================================================
// d(vol)w / dJ
//================================================================================
const double  SimoPisterWEnergy::dvolwOdJ(const double &J_in )
{
   double dcolwOverdJ = K * (-2.0 / J_in + 2.0 * J_in) * 0.25;
   return dcolwOverdJ;
}

//================================================================================
// d2(vol)w / dJ2
//================================================================================
const double  SimoPisterWEnergy::d2volwOdJ2(const double &J_in )
{
   double d2colwOverdJ2 = K * (2.0 / J_in / J_in + 2.0) * 0.25;
   return d2colwOverdJ2;
}


#endif

