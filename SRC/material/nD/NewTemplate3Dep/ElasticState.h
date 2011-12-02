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

#ifndef ElasticState_H
#define ElasticState_H

#include "iostream.h"
#include <BJtensor.h>
#include <stresst.h>
#include <straint.h>
#include "MaterialParameter.h"

class ElasticState
{  
  public:
    
    ElasticState(const stresstensor &initialStress, const straintensor &initialStrain);
    ElasticState(const stresstensor &initialStress);
    ElasticState();
    
    virtual ~ElasticState() {};
    virtual ElasticState* newObj() = 0;

    virtual stresstensor getStress() const;
    virtual straintensor getStrain() const;
    
    virtual const BJtensor &getElasticStiffness (const MaterialParameter &MatPar_in) const = 0;
    
    virtual int setStress(const stresstensor &Stre_in);
    virtual int setStrain(const straintensor &Stra_in);
  
  protected:                 
    
    stresstensor Stress;
    straintensor Strain; 
    
    static BJtensor ElasticStiffness;
    static const stresstensor zerostress;
    static const straintensor zerostrain;

};


#endif

