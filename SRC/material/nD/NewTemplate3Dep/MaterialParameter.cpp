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

// This is the container to store material constants:
// Material_Parameter:      to store fixed scalar material constants;
//                          Note: the first one should be the density, input 0.0 if not involved;
// Num_Material_Parameter:  the number of Material_Parameter, should not less than 2, 
//                          given the 1st is the density and the 2nd is void ratio;
// Internal_Scalar:         to store initial evolution scalar variables first time, 
//                          and scalar variables thereafter;
// Num_Internal_Scalar:     the number of Internal_Scalar;
// Internal_Tensor:         to store initial evolution tensor variables first time (usually zero tensor), 
//                          and tensor variables thereafter;
// Num_Internal_Tensor:     the number of Internal_Tensor;


#ifndef MaterialParameter_CPP
#define MaterialParameter_CPP

#include "MaterialParameter.h"

// Constructor 1
MaterialParameter::MaterialParameter(const double *Material_Parameter_in, 
                                    int Num_Material_Parameter_in, 
                                    const double *Internal_Scalar_in, 
                                    int Num_Internal_Scalar_in, 
                                    const stresstensor *Internal_Tensor_in, 
                                    int Num_Internal_Tensor_in)
{
    // Material Constants
    if (Num_Material_Parameter_in > 0) {
        Material_Parameter = new double [Num_Material_Parameter = Num_Material_Parameter_in];
        if (!Material_Parameter) {
            cout << "MaterialParameter::Insufficient memory! " << endl;
            exit (1);
        }        
        for (int i = 0; i < Num_Material_Parameter; i++)
            Material_Parameter[i] = Material_Parameter_in[i];
    }
    else {
      Num_Material_Parameter = 0;
      Material_Parameter = NULL;
    }

    // Scalar (Isotropic) Internal Variables
    if (Num_Internal_Scalar_in > 0) {
        Internal_Scalar = new double [Num_Internal_Scalar = Num_Internal_Scalar_in];
        if (!Internal_Scalar) {
            cout << "MaterialParameter::Insufficient memory! " << endl;
            exit (1);
        }        
        for (int i = 0; i < Num_Internal_Scalar; i++)
            Internal_Scalar[i] = Internal_Scalar_in[i];
    }
    else {
      Num_Internal_Scalar = 0;
      Internal_Scalar = NULL;
    }

    // Tensor (Kinematic) Internal Variables
    if (Num_Internal_Tensor_in > 0) {
        Internal_Tensor = new stresstensor[Num_Internal_Tensor = Num_Internal_Tensor_in];
        if (!Internal_Tensor) {
            cout << "MaterialParameter::Insufficient memory! " << endl;
            exit (1);
        }        
        for (int i = 0; i < Num_Internal_Tensor; i++)
            Internal_Tensor[i] = Internal_Tensor_in[i];
    }
    else {
      Num_Internal_Tensor = 0;
      Internal_Tensor = NULL;
    }
}

// Constructor 2
MaterialParameter::MaterialParameter(const double *Material_Parameter_in, 
                                     int Num_Material_Parameter_in, 
                                     const stresstensor *Internal_Tensor_in, 
                                     int Num_Internal_Tensor_in)
: Internal_Scalar(NULL), Num_Internal_Scalar(0)
{
    // Material Constants
    if (Num_Material_Parameter_in > 0) {
        Material_Parameter = new double [Num_Material_Parameter = Num_Material_Parameter_in];
        if (!Material_Parameter) {
            cout << "MaterialParameter::Insufficient memory! " << endl;
            exit (1);
        }        
        for (int i = 0; i < Num_Material_Parameter; i++)
            Material_Parameter[i] = Material_Parameter_in[i];
    }
    else {
      Num_Material_Parameter = 0;
      Material_Parameter = NULL;
    }

    // Tensor (Kinematic) Internal Variables
    if (Num_Internal_Tensor_in > 0) {
        Internal_Tensor = new stresstensor[Num_Internal_Tensor = Num_Internal_Tensor_in];
        if (!Internal_Tensor) {
            cout << "MaterialParameter::Insufficient memory! " << endl;
            exit (1);
        }        
        for (int i = 0; i < Num_Internal_Tensor; i++)
            Internal_Tensor[i] = Internal_Tensor_in[i];
    }
    else {
      Num_Internal_Tensor = 0;
      Internal_Tensor = NULL;
    }
}

// Constructor 3
MaterialParameter::MaterialParameter( )
: Material_Parameter(NULL), Num_Material_Parameter(0), 
  Internal_Scalar(NULL), Num_Internal_Scalar(0), 
  Internal_Tensor(NULL), Num_Internal_Tensor(0)
{

}

// Destructor
MaterialParameter::~MaterialParameter( )
{
	if (Material_Parameter)
		delete [] Material_Parameter;

	if (Internal_Scalar)
		delete [] Internal_Scalar;

	if (Internal_Tensor)
		delete [] Internal_Tensor;
}

// Copy constructor
MaterialParameter::MaterialParameter(const MaterialParameter &refer_MaterialParameter )
{
    // Material Constants
    if (refer_MaterialParameter.getNum_Material_Parameter() > 0) {
        Material_Parameter = new double [Num_Material_Parameter = refer_MaterialParameter.getNum_Material_Parameter()];
        if (!Material_Parameter) {
            cout << "MaterialParameter::Insufficient memory! " << endl;
            exit (1);
        }        
        for (int i = 0; i < Num_Material_Parameter; i++)
            Material_Parameter[i] = refer_MaterialParameter.getMaterial_Parameter(i);
    }
    else {
      Num_Material_Parameter = 0;
      Material_Parameter = NULL;
    }

    // Scalar (Isotropic) Internal Variables
    if (refer_MaterialParameter.getNum_Internal_Scalar() > 0) {
        Internal_Scalar = new double [Num_Internal_Scalar = refer_MaterialParameter.getNum_Internal_Scalar()];
        if (!Internal_Scalar) {
            cout << "MaterialParameter::Insufficient memory! " << endl;
            exit (1);
        }        
        for (int i = 0; i < Num_Internal_Scalar; i++)
            Internal_Scalar[i] = refer_MaterialParameter.getInternal_Scalar(i);
    }
    else {
      Num_Internal_Scalar = 0;
      Internal_Scalar = NULL;
    }

    // Tensor (Kinematic) Internal Variables
    if (refer_MaterialParameter.getNum_Internal_Tensor() > 0) {
        Internal_Tensor = new stresstensor[Num_Internal_Tensor = refer_MaterialParameter.getNum_Internal_Tensor()];
        if (!Internal_Tensor) {
            cout << "MaterialParameter::Insufficient memory! " << endl;
            exit (1);
        }        
        for (int i = 0; i < Num_Internal_Tensor; i++)
            Internal_Tensor[i] = refer_MaterialParameter.getInternal_Tensor(i);
    }
    else {
      Num_Internal_Tensor = 0;
      Internal_Tensor = NULL;
    }

}

// Create a new class pointer
MaterialParameter* MaterialParameter::newObj() {
	MaterialParameter *ptr_MaterialParameter = new  MaterialParameter ( this->Material_Parameter,
                                                                        this->Num_Material_Parameter,
                                                                        this->Internal_Scalar,
                                                                        this->Num_Internal_Scalar,
                                                                        this->Internal_Tensor,
                                                                        this->Num_Internal_Tensor );
    return ptr_MaterialParameter;	
}

// getNum_Material_Parameter
int MaterialParameter::getNum_Material_Parameter() const
{
	return Num_Material_Parameter;
}

// getNum_Internal_Scalar
int MaterialParameter::getNum_Internal_Scalar() const
{
	return Num_Internal_Scalar;
}

// getNum_Internal_Tensor
int MaterialParameter::getNum_Internal_Tensor() const
{
	return Num_Internal_Tensor;
}

// getMaterial_Parameter
double MaterialParameter::getMaterial_Parameter(int which) const
{
	if (which < 0 || which > Num_Material_Parameter) {
		cout << "Error! MaterialParameter::getMaterial_Parameter - Invalid index of material constants. " << endl;
		exit(1);
	}	

	return Material_Parameter[which];
}

// getInternal_Scalar
double MaterialParameter::getInternal_Scalar(int which) const
{
	if (which < 0 || which > Num_Internal_Scalar) {
		cout << "Error! MaterialParameter::getInternal_Scalar - Invalid index of internal scalars. " << endl;
		exit(1);
	}	

	return Internal_Scalar[which];
}

// getInternal_Tensor
const stresstensor& MaterialParameter::getInternal_Tensor(int which) const
{
	if (which < 0 || which > Num_Internal_Tensor) {
		cout << "Error! MaterialParameter::getInternal_Tensor - Invalid index of internal tensors. " << endl;
		exit (1);
	}	

	return Internal_Tensor[which];
}

// setMaterial_Parameter
int MaterialParameter::setMaterial_Parameter(int which, double newMaterial_Parameter) 
{
	if (which < 0 || which > Num_Material_Parameter) {
		cout << "Error! MaterialParameter::setMaterial_Parameter - Invalid index of material constants. " << endl;
		return (1);
	}	
	
	Material_Parameter[which] = newMaterial_Parameter;
	
	return 0;
}

// setInternal_Scalar
int MaterialParameter::setInternal_Scalar(int which, double newInternal_Scalar) 
{
	if (which < 0 || which > Num_Internal_Scalar) {
		cout << "Error! MaterialParameter::setInternal_Scalar - Invalid index of internal scalars. " << endl;
		exit (1);
	}	
	
	Internal_Scalar[which] = newInternal_Scalar;
	
	return 0;
}

// setInternal_Tensor
int MaterialParameter::setInternal_Tensor(int which, const stresstensor &newInternal_Tensor) 
{
	if (which < 0 || which > Num_Internal_Tensor) {
		cout << "Error! MaterialParameter::setInternal_Tensor - Invalid index of internal tensors. " << endl;
		exit (1);
	}	
		
	Internal_Tensor[which] = newInternal_Tensor;
				
	return 0;
}


#endif

