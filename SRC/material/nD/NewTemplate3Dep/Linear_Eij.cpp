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

#ifndef Linear_Eij_CPP
#define Linear_Eij_CPP

#include "Linear_Eij.h"
#include <Channel.h>
#include <ID.h>

Linear_Eij::Linear_Eij(int LinearFactor_index_in)
  :TensorEvolution(TENSOR_EVOLUTION_TAGS_Linear_Eij),
   LinearFactor_index(LinearFactor_index_in)
{

}

TensorEvolution* Linear_Eij::newObj()
{
   TensorEvolution* nObj = new Linear_Eij(this->LinearFactor_index);

   return nObj;
}

const straintensor& Linear_Eij::Hij(const straintensor& plastic_flow, const stresstensor& Stre, 
                            const straintensor& Stra, const MaterialParameter& material_parameter)
{
    TensorEvolution::TensorEvolutionHij = plastic_flow * getLinearFactor(material_parameter); 
    return TensorEvolution::TensorEvolutionHij;
}

double Linear_Eij::getLinearFactor(const MaterialParameter& material_parameter) const
{
    if ( LinearFactor_index <= material_parameter.getNum_Material_Parameter() && LinearFactor_index > 0)
        return material_parameter.getMaterial_Parameter(LinearFactor_index -1);
    else {
        opserr << "Linear_Eij: Invalid Input. " << endln;
        exit (1);
    }  
}

int 
Linear_Eij::sendSelf(int commitTag, Channel &theChannel)
{
  static ID iData(1);
  iData(0) = LinearFactor_index;

  int dbTag = this->getDbTag();

  if (theChannel.sendID(dbTag, commitTag, iData) < 0) {
    opserr << "Linear_Eij::sendSelf() - failed to send data\n";
    return -1;
  }

  return 0;
}

int 
Linear_Eij::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID iData(1);
  int dbTag = this->getDbTag();

  if (theChannel.recvID(dbTag, commitTag, iData) < 0) {
    opserr << "Linear_Eij::recvSelf() - failed to recv data\n";
    return -1;
  }

  LinearFactor_index = iData(0);

  return 0;
}

#endif

