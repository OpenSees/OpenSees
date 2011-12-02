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

#ifndef NewTemplate3Dep_H
#define NewTemplate3Dep_H

#include <stresst.h>
#include <straint.h>
#include <BJtensor.h>
#include <iostream.h>

#include <NDMaterial.h>

#include "MaterialParameter.h"
#include "ElasticState.h"
#include "YieldFunction.h"
#include "PlasticFlow.h"
#include "ScalarEvolution.h"
#include "TensorEvolution.h"

#include <Channel.h>
#include <G3Globals.h>

class NewTemplate3Dep : public NDMaterial
{
public:
    NewTemplate3Dep( int tag,
                  MaterialParameter *pointer_material_parameter_in,
                  ElasticState      *pointer_elastic_state_in,
                  YieldFunction     *pointer_yield_function_in,        
                  PlasticFlow       *pointer_plastic_flow_in,
                  ScalarEvolution   **pointer_scalar_evolution_in = NULL,
                  TensorEvolution   **pointer_tensor_evolution_in = NULL,
                  int caseIndex_in = 0);
    
     NewTemplate3Dep( int tag,
                  MaterialParameter *pointer_material_parameter_in,
                  ElasticState      *pointer_elastic_state_in,
                  YieldFunction     *pointer_yield_function_in,        
                  PlasticFlow       *pointer_plastic_flow_in,
                  TensorEvolution   **pointer_tensor_evolution_in,
                  int caseIndex_in = 0);
    
    NewTemplate3Dep(void);
    
    ~NewTemplate3Dep(void);

    // methods to set and retrieve state using the Tensor class    
    int setTrialStrain(const Tensor& v);
    int setTrialStrain(const Tensor& v, const Tensor& r);    
    int setTrialStrainIncr(const Tensor& v);
    int setTrialStrainIncr(const Tensor& v, const Tensor& r);
    
    double getRho();

    const BJtensor& getTangentTensor(void);
    const stresstensor& getStressTensor(void);
    const straintensor& getStrainTensor(void);
    const straintensor& getPlasticStrainTensor(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    
    NDMaterial* getCopy(void);
    NDMaterial* getCopy(const char *code);

    const char *getType(void) const;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

    void Print(OPS_Stream& s, int flag =0);
    //friend OPS_Stream& operator<< (OPS_Stream& os, const NewTemplate3Dep& MP);                                              
  
private:

    int ForwardEuler(const straintensor& strain_incr);
    int SemiImplicit(const straintensor& strain_incr);
    int BackwardEuler(const straintensor& strain_incr); 

    int PredictorEPState(const straintensor& strain_incr);
    stresstensor yield_surface_cross(const stresstensor& start_stress, 
                                     const stresstensor& end_stress,
				     double a);
    double zbrentstress(const stresstensor& start_stress,
                        const stresstensor& end_stress,
                        double x1, double x2, double tol) const;
    double func( const stresstensor& start_stress,
                 const stresstensor& end_stress,
                 const MaterialParameter& pointer_material_parameter, 
                 double alfa ) const;
 
private:
   
    straintensor TrialStrain;
    stresstensor TrialStress;
    straintensor TrialPlastic_Strain;
    
    stresstensor CommitStress;
    straintensor CommitStrain; 
    straintensor CommitPlastic_Strain;
    
    BJtensor Stiffness;
                
    MaterialParameter *pointer_material_parameter;
    ElasticState      *pointer_elastic_state;           
    YieldFunction     *pointer_yield_function;
    PlasticFlow       *pointer_plastic_flow;
    ScalarEvolution   **pointer_scalar_evolution;
    TensorEvolution   **pointer_tensor_evolution;
    
    int caseIndex;
};

#endif
