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


#ifndef W_H
#define W_H

#include <Vector.h>
#include <Tensor.h>

class WEnergy
{
  public:
    WEnergy();
    virtual ~WEnergy();
    virtual WEnergy*newObj( ) =0;

    virtual const  double wE(const double &, const Vector &) ;
    virtual const Vector disowOdlambda(const Vector &) ;
    virtual const Vector d2isowOdlambda2(const Vector &) ;
    virtual const Tensor d2isowOdlambda1dlambda2(const Vector &) ;
    virtual const double dvolwOdJ(const double &) ;
    virtual const double d2volwOdJ2(const double &) ;

};

#endif

