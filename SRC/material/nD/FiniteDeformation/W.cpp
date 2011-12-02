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
//# UPDATE HISTORY:
//#
//#
//===============================================================================


#ifndef W_cpp
#define W_cpp


#include "W.h"

WEnergy::WEnergy()
{

}

WEnergy::~WEnergy()
{

}

const double WEnergy::wE(const double &J_in, const Vector &lambda_wave_in)
{
  return 0.0;
}


const Vector WEnergy::disowOdlambda( const Vector &lambda_wave_in)
{
  return Vector(3);
}

const Vector WEnergy::d2isowOdlambda2( const Vector &lambda_wave_in)
{
  return Vector(3);
}

const Tensor WEnergy::d2isowOdlambda1dlambda2(const Vector &lambda_wave_in)
{
  Tensor zerotensor(2,def_dim_2,0.0);
  return zerotensor;
}

const double WEnergy::dvolwOdJ(const double &J_in)
{
  return 0.0;
}

const double WEnergy::d2volwOdJ2(const double &J_in)
{
  return 0.0;;

}


#endif

