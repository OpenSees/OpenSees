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
//# UPDATE HISTORY:    28May2004
//#
//#
//===============================================================================

#ifndef OgdenWEnergy_H
#define OgdenWEnergy_H

#include <Vector.h>
#include <Tensor.h>
#include <OPS_Globals.h>
#include <W.h>

#include <math.h>


class OgdenWEnergy : public WEnergy
{
  private:
    int N_Ogden;
    double *cr_Ogden;
    double *mur_Ogden;
  public:
    OgdenWEnergy(int , double * , double * );
    OgdenWEnergy( );
    ~OgdenWEnergy( );
    WEnergy *newObj( );

    const double  wE(const double &, const Vector &) ;
    const Vector  disowOdlambda(const Vector & )  ;
    const Vector  d2isowOdlambda2(const Vector & )  ;
//    const tensor  d2isowOdlambda1dlambda2( const Vector &)  ;
//    const double  dvolwOdJ( const double &)  ;
//    const double  d2volwOdJ2( const double &) ;


};

#endif

