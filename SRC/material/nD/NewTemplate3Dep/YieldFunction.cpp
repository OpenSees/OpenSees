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

#ifndef YieldFunction_CPP
#define YieldFunction_CPP

#include "YieldFunction.h"

double YieldFunction::InScalarDerivative(const stresstensor& Stre, 
                                         const MaterialParameter &MaterialParameter_in, 
                                         int which) const
{
    return 0.0;
}

const stresstensor& YieldFunction::InTensorDerivative(const stresstensor& Stre, 
                                                      const MaterialParameter &MaterialParameter_in, 
                                                      int which) const
{
    // Just for compiling
    static stresstensor ZStre;
    return ZStre;
}

//int YieldFunction::getTensionOrCompressionType() const
//{
//    // if No Tensile allowed, return 1;
//    // By default, return 0;
//    return 0;
//}

#endif

