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

#ifndef TensorEvolution_CPP
#define TensorEvolution_CPP

#include "TensorEvolution.h"

straintensor TensorEvolution::TensorEvolutionHij;

TensorEvolution::TensorEvolution()
{

}

//TensorEvolution* TensorEvolution::newObj()
//{
//	TensorEvolution* nObj = new TensorEvolution();
//	
//	return nObj;
//}

const straintensor& TensorEvolution::Hij(const straintensor& plastic_flow, const stresstensor& Stre, 
                                         const straintensor& Stra, const MaterialParameter& material_parameter)
{   
    return TensorEvolution::TensorEvolutionHij;
}


#endif

