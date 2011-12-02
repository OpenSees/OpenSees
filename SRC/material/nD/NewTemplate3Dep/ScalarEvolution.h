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

#ifndef ScalarEvolution_H
#define ScalarEvolution_H 

#include <straint.h>
#include <stresst.h>
#include "MaterialParameter.h"
#include "ElasticState.h"
#include <MovableObject.h>

class ScalarEvolution : public MovableObject
{
  public:
    
    ScalarEvolution(int classTag);
    virtual ~ScalarEvolution() {};

    virtual ScalarEvolution *newObj() = 0;

    virtual double H(const straintensor& plastic_flow, const stresstensor& Stre, 
                     const straintensor& Stra, const MaterialParameter& material_parameter);
                     
};


#endif

