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

#ifndef MaterialParameter_H
#define MaterialParameter_H

#include <stresst.h>
#include <iostream.h>
		
class MaterialParameter
{  
  public:

    MaterialParameter(const double *Material_Parameter_in = NULL, 
                      int Num_Material_Parameter_in = 0, 
		      const double *Internal_Scalar_in = NULL, 
		      int Num_Internal_Scalar_in = 0, 
		      const stresstensor *Internal_Tensor_in = NULL, 
		      int Num_Internal_Tensor_in = 0);
    MaterialParameter(const double *Material_Parameter_in, 
                      int Num_Material_Parameter_in, 
		      const stresstensor *Internal_Tensor_in, 
		      int Num_Internal_Tensor_in);
    MaterialParameter( );
    ~MaterialParameter( );
    MaterialParameter(const MaterialParameter &refer_MaterialParameter);
    MaterialParameter* newObj();

    int getNum_Material_Parameter() const;
    int getNum_Internal_Scalar() const;
    int getNum_Internal_Tensor() const;
    double getMaterial_Parameter(int which) const;
    double getInternal_Scalar(int which) const;
    const stresstensor& getInternal_Tensor(int which) const;
    
    int setMaterial_Parameter(int which, double newMaterial_Parameter);
    int setInternal_Scalar(int which, double newInternal_Scalar);
    int setInternal_Tensor(int which, const stresstensor &newInternal_Tensor);    
  
  private:
                 
    double *Material_Parameter; 
    int Num_Material_Parameter;
    double *Internal_Scalar;    
    int Num_Internal_Scalar;
    stresstensor *Internal_Tensor;  
    int Num_Internal_Tensor;
};


#endif

