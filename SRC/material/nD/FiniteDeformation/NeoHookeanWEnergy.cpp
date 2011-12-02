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

#ifndef NeoHookeanWEnergy_CPP
#define NeoHookeanWEnergy_CPP

#include <NeoHookeanWEnergy.h>

//================================================================================
// Normal constructor
//================================================================================
NeoHookeanWEnergy::NeoHookeanWEnergy(double K_in, double G_in )
 :K(K_in), G(G_in)
{

}

NeoHookeanWEnergy::NeoHookeanWEnergy(  )
 :K(0.0), G(0.0)
{

}

//================================================================================
// Normal destructor
//================================================================================
NeoHookeanWEnergy::~NeoHookeanWEnergy( )
{

}

//================================================================================
//create a clone of itself
//================================================================================
WEnergy * NeoHookeanWEnergy::newObj( )
  {
    WEnergy  *new_WEnergy = new NeoHookeanWEnergy(K,  G);
    return new_WEnergy;
  }


//================================================================================
// w
//================================================================================
const double  NeoHookeanWEnergy::wE(const double &J_in, const Vector &lambda_wave_in )
  {
    double w_iso = 0.5 * G * (lambda_wave_in(0) * lambda_wave_in(0)
                            + lambda_wave_in(1) * lambda_wave_in(1)
                            + lambda_wave_in(2) * lambda_wave_in(2) );

    double w_vol = 0.5 * K * (J_in-1.0) * (J_in-1.0); //version I
//    double w_vol = 0.25 * K * ( J_in*J_in -1.0 -2.0*log(J_in) );  //version II
//    double w_vol = 0.0; //incompressible material
    double w_total = w_iso + w_vol;
    return w_total;
  }

//================================================================================
// d(iso)w / d(lambda)
//================================================================================
const Vector  NeoHookeanWEnergy::disowOdlambda(const Vector &lambda_wave_in )
  {
    Vector disowOverdlambda(3);
    disowOverdlambda(0) = G * lambda_wave_in(0) ;
    disowOverdlambda(1) = G * lambda_wave_in(1) ;
    disowOverdlambda(2) = G * lambda_wave_in(2) ;

    return disowOverdlambda;
  }

//================================================================================
// d2(iso)w / d(lambda)2
//================================================================================
const Vector  NeoHookeanWEnergy::d2isowOdlambda2(const Vector &lambda_wave_in )
  {
    Vector d2isowOverdlambda2(3);
    d2isowOverdlambda2(0) = G;
    d2isowOverdlambda2(1) = G;
    d2isowOverdlambda2(2) = G;
    return d2isowOverdlambda2;
  }


//================================================================================
// d(vol)w / dJ
//================================================================================
const double  NeoHookeanWEnergy::dvolwOdJ(const double &J_in )
{
  double temp1 = K * (J_in - 1.0);  // Version I
//  double temp1 = 0.0;
  return temp1;
}

//================================================================================
// d2(vol)w / dJ2
//================================================================================
const double  NeoHookeanWEnergy::d2volwOdJ2(const double &J_in )
{
  double temp2 = K ;  // version I
//  double temp2 = 0.0;
  return temp2;
}


#endif

