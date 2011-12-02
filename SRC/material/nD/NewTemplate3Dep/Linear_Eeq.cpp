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

#ifndef Linear_Eeq_CPP
#define Linear_Eeq_CPP

#include "Linear_Eeq.h"
#include <Channel.h>
#include <ID.h>

Linear_Eeq::Linear_Eeq(int LinearFactor_index_in)
  : ScalarEvolution(SCALAR_EVOLUTION_TAGS_Linear_Eeq), 
    LinearFactor_index(LinearFactor_index_in)
{

}

ScalarEvolution* Linear_Eeq::newObj()
{
    ScalarEvolution* nObj = new Linear_Eeq(this->LinearFactor_index);
    
    return nObj;
}

double Linear_Eeq::H(const straintensor& plastic_flow, const stresstensor& Stre, 
                     const straintensor& Stra, const MaterialParameter& material_parameter)
{
    return plastic_flow.equivalent() * getLinearFactor(material_parameter);
}

double Linear_Eeq::getLinearFactor(const MaterialParameter& material_parameter) const
{
    if ( LinearFactor_index <= material_parameter.getNum_Material_Parameter() && LinearFactor_index > 0)
        return material_parameter.getMaterial_Parameter(LinearFactor_index-1);
    else {
        opserr << "Linear_Eeq: Invalid Input. " << endln;
        exit (1);
    }
} 

int 
Linear_Eeq::sendSelf(int commitTag, Channel &theChannel)
{
  static ID iData(1);
  iData(0) = LinearFactor_index;
  int dbTag = this->getDbTag();

  if (theChannel.sendID(dbTag, commitTag, iData) < 0) {
    opserr << "Linear_Eeq::sendSelf() - failed to send data\n";
    return -1;
  }

  return 0;
}
int 
Linear_Eeq::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID iData(1);
  int dbTag = this->getDbTag();

  if (theChannel.recvID(dbTag, commitTag, iData) < 0) {
    opserr << "Linear_Eeq::recvSelf() - failed to recv data\n";
    return -1;
  }

  LinearFactor_index = iData(0);

  return 0;
}

#endif

