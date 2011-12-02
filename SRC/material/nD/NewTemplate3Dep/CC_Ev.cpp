///////////////////////////////////////////////////////////////////////////////
//   COPYLEFT (C): Woody's viral GPL-like license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              
// CLASS:             
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           
//
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         
// DESIGNER:          Zhao Cheng, Boris Jeremic
// PROGRAMMER:        Zhao Cheng, 
// DATE:              Fall 2005
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef CC_Ev_CPP
#define CC_Ev_CPP

#include "CC_Ev.h"

CC_Ev::CC_Ev(int lambda_index_in, 
             int kappa_index_in,
             int e0_index_in,
             int p0_index_in)
: lambda_index(lambda_index_in), 
  kappa_index(kappa_index_in),
  e0_index(e0_index_in),
  p0_index(p0_index_in)
{

}

ScalarEvolution* CC_Ev::newObj()
{
    ScalarEvolution* nObj = new CC_Ev(this->lambda_index,
                                      this->kappa_index,
                                      this->e0_index,
                                      this->p0_index);
    return nObj;
}

double CC_Ev::H(const straintensor& plastic_flow, const stresstensor& Stre, 
                const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double d_Ev = - plastic_flow.Iinvariant1(); // note "minus"
    double p0 = getp0(material_parameter);
    double lambda = getlambda(material_parameter);
    double kappa = getkappa(material_parameter);
    double e0 = gete0(material_parameter);
    
    double e = e0 + (1.0 + e0) *Stra.Iinvariant1();

    return (1.0 + e) * p0 * d_Ev / (lambda - kappa);
}

// Get lambda
double CC_Ev::getlambda(const MaterialParameter& material_parameter) const
{
    if ( lambda_index <= material_parameter.getNum_Material_Parameter() && lambda_index > 0)
        return material_parameter.getMaterial_Parameter(lambda_index-1);
    else {
        cout << "CC_Ev: Invalid Input of " << lambda_index << endl;
        exit (1);
    }
}

// Get kappa
double CC_Ev::getkappa(const MaterialParameter& material_parameter) const
{
    if ( kappa_index <= material_parameter.getNum_Material_Parameter() && kappa_index > 0)
        return material_parameter.getMaterial_Parameter(kappa_index-1);
    else {
        cout << "CC_Ev: Invalid Input of " << kappa_index << endl;
        exit (1);
    }
}

// Get e0
double CC_Ev::gete0(const MaterialParameter& material_parameter) const
{
    if ( e0_index <= material_parameter.getNum_Material_Parameter() && e0_index > 0)
        return material_parameter.getMaterial_Parameter(e0_index-1);
    else {
        cout << "CC_Ev: Invalid Input of " << e0_index << endl;
        exit (1);
    }
}

// Get p0
double CC_Ev::getp0(const MaterialParameter& material_parameter) const
{
    if ( p0_index <= material_parameter.getNum_Internal_Scalar() && p0_index > 0)
        return material_parameter.getInternal_Scalar(p0_index-1);
    else {
        cout << "CC_Ev: Invalid Input of " << p0_index << endl;
        exit (1);
    }
}

#endif

