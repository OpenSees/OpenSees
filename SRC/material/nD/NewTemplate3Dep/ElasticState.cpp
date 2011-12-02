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

#ifndef ElasticState_CPP
#define ElasticState_CPP

#include "ElasticState.h"

BJtensor ElasticState::ElasticStiffness(4, def_dim_4, 0.0);
const stresstensor ElasticState::zerostress(0.0);
const straintensor ElasticState::zerostrain(0.0);

////////////////////////////////////////////////////////////////
ElasticState::ElasticState(const stresstensor &initialStress, const straintensor &initialStrain)
: Stress(initialStress), Strain(initialStrain)
{

}

////////////////////////////////////////////////////////////////
ElasticState::ElasticState(const stresstensor &initialStress)
: Stress(initialStress)
{
    straintensor ZeroStra;
    Strain = ZeroStra;
}

////////////////////////////////////////////////////////////////
ElasticState::ElasticState()
{
    stresstensor ZeroStre;
    Stress = ZeroStre;
    
    straintensor ZeroStra;
    Strain = ZeroStra;
}

////////////////////////////////////////////////////////////////
stresstensor ElasticState::getStress() const 
{ 
    return Stress;
}

////////////////////////////////////////////////////////////////
straintensor ElasticState::getStrain() const 
{ 
    return Strain;
}

/////////////////////////////////////////////////////////////////
int ElasticState::setStress(const stresstensor &Stre_in) 
{
    Stress = Stre_in;
    
    return 0;
}

/////////////////////////////////////////////////////////////////
int ElasticState::setStrain(const straintensor &Stra_in) 
{
    Strain = Stra_in;
    
    return 0;
}

#endif
