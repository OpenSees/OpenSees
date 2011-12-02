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

Linear_Eeq::Linear_Eeq(int LinearFactor_index_in)
: LinearFactor_index(LinearFactor_index_in)
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
        cout << "Linear_Eeq: Invalid Input. " << endl;
        exit (1);
    }
} 


#endif

