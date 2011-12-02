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

// Ref: Dafalias and Manzari 2004: J. Eng. Mech. 130(6), pp 622-634
// Parameters:
//  1- m:         parameter in the yield function;
//  2- c_z:       parameter
//  3- z_max      parameter
//  4- alpha:     "back-stress" tensor in yield function; (the 1st tensorial internal variable);
//  5- z:         fabric dilation internal tensor (the 2nd tensorial internal variable); 

#ifndef DM04_z_Eij_CPP
#define DM04_z_Eij_CPP

#include "DM04_z_Eij.h"
#include <Channel.h>
#include <ID.h>

stresstensor DM04_z_Eij::DM04_z_t;

DM04_z_Eij::DM04_z_Eij(int m_index_in,
                       int c_z_index_in,
                       int z_max_index_in,
                       int alpha_index_in,
                       int z_index_in)
  :TensorEvolution(TENSOR_EVOLUTION_TAGS_DM04_z_Eij),
   m_index(m_index_in), 
  c_z_index(c_z_index_in), 
  z_max_index(z_max_index_in),
  alpha_index(alpha_index_in), 
  z_index(z_index_in)
{

}

TensorEvolution* DM04_z_Eij::newObj()
{
    TensorEvolution* nObj = new DM04_z_Eij(this->m_index,
                                           this->c_z_index,
                                           this->z_max_index,
                                           this->alpha_index,
                                           this->z_index);

    return nObj;
}

const straintensor& DM04_z_Eij::Hij(const straintensor& plastic_flow, const stresstensor& Stre, 
                                    const straintensor& Stra, const MaterialParameter& material_parameter)
{
    double m = getm(material_parameter);
    double c_z = getc_z(material_parameter);
    double z_max = getz_max(material_parameter);        
    stresstensor alpha = getalpha(material_parameter);
    stresstensor z = getz(material_parameter);
    
    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();

    stresstensor n;
    
    stresstensor s_bar = Stre.deviator() - (alpha *p);
    double _s_bar_ = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    if (p > 0.0 && _s_bar_ > 0.0)
        n = s_bar * (1.0/_s_bar_);
    
    // here d_Ev has different sign assumption from the ref.
    // hence no "negative" sign for d_Ev in the following line
    double d_Ev = plastic_flow.Iinvariant1();
    if (d_Ev < 0.0) 
      d_Ev = 0.0;  
   
    TensorEvolution::TensorEvolutionHij = ((n *z_max) +z) *(-c_z*d_Ev); 
    
    return TensorEvolution::TensorEvolutionHij;
}

// to get m
//================================================================================
double DM04_z_Eij::getm(const MaterialParameter& material_parameter) const
{
    if ( m_index <= material_parameter.getNum_Material_Parameter() && m_index > 0)
        return material_parameter.getMaterial_Parameter(m_index-1);
    else {
        opserr << "DM04_alpha: Invalid Input. " << endln;
        exit (1);
    }
}

// to get c_z
//================================================================================
double DM04_z_Eij::getc_z(const MaterialParameter& material_parameter) const
{
    if ( c_z_index <= material_parameter.getNum_Material_Parameter() && c_z_index > 0)
        return material_parameter.getMaterial_Parameter(c_z_index-1);
    else {
        opserr << "DM04_alpha: Invalid Input. " << endln;
        exit (1);
    }
}

// to get c
//================================================================================
double DM04_z_Eij::getz_max(const MaterialParameter& material_parameter) const
{
    if ( z_max_index <= material_parameter.getNum_Material_Parameter() && z_max_index > 0)
        return material_parameter.getMaterial_Parameter(z_max_index-1);
    else {
        opserr << "DM04_alpha: Invalid Input. " << endln;
        exit (1);
    }
}

// to get alpha
//================================================================================
const stresstensor& DM04_z_Eij::getalpha(const MaterialParameter& material_parameter) const
{
    if ( alpha_index <= material_parameter.getNum_Internal_Tensor() && alpha_index > 0) {
        DM04_z_Eij::DM04_z_t = material_parameter.getInternal_Tensor(alpha_index-1);
        return DM04_z_Eij::DM04_z_t;
    }
    else {
        opserr << "DM04_z: Invalid Input. " << endln;
        exit (1);
    }
}

// to get z
//================================================================================
const stresstensor& DM04_z_Eij::getz(const MaterialParameter& material_parameter) const
{
    if ( z_index <= material_parameter.getNum_Internal_Tensor() && z_index > 0) {
        DM04_z_Eij::DM04_z_t = material_parameter.getInternal_Tensor(z_index-1);
        return DM04_z_Eij::DM04_z_t;
    }
    else {
        opserr << "DM04_z: Invalid Input. " << endln;
        exit (1);
    }
}

int
DM04_z_Eij::sendSelf(int commitTag, Channel &theChannel)
{
  static ID iData(5);

  iData(0) = m_index;
  iData(1) = c_z_index;
  iData(2) = z_max_index;
  iData(3) = alpha_index;
  iData(4) = z_index;

  int dbTag = this->getDbTag();

  if (theChannel.sendID(dbTag, commitTag, iData) < 0) {
    opserr << "DM04_z_Eij::sendSelf() - failed to send data\n";
    return -1;
  }

  return 0;
}

int 
DM04_z_Eij::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID iData(5);
  int dbTag = this->getDbTag();

  if (theChannel.recvID(dbTag, commitTag, iData) < 0) {
    opserr << "DM04_z_Eij::recvSelf() - failed to recv data\n";
    return -1;
  }

  m_index = iData(0);
  c_z_index = iData(1);
  z_max_index = iData(2);
  alpha_index = iData(3);
  z_index = iData(4);


  return 0;
}

#endif

