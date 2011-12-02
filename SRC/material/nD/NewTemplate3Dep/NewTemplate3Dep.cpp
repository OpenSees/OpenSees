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
// UPDATE HISTORY:    06/2006, add functions for matrix based elements, CZ
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef NewTemplate3Dep_CPP
#define NewTemplate3Dep_CPP

#include "NewTemplate3Dep.h"

const  straintensor NewTemplate3Dep::ZeroStrain;
const  stresstensor NewTemplate3Dep::ZeroStress;
const  BJtensor NewTemplate3Dep::ZeroI4(4, def_dim_4, 0.0);
const int NewTemplate3Dep::ISMAX = 30;
const int NewTemplate3Dep::ITMAX = 30;
const double NewTemplate3Dep::TOL = 1.0e-7;
const double NewTemplate3Dep::FTOL = 1.0e-8;

// For Matrix based elements
Matrix NewTemplate3Dep::D(6,6);
Vector NewTemplate3Dep::sigma(6);
Vector NewTemplate3Dep::epsilon(6);

#include "NewTemplate3Dep.h"

// Constructor
//================================================================================
NewTemplate3Dep::NewTemplate3Dep( int tag,
                                  MaterialParameter *pointer_material_parameter_in,
                                  ElasticState      *pointer_elastic_state_in,
                                  YieldFunction     *pointer_yield_function_in ,        
                                  PlasticFlow       *pointer_plastic_flow_in,
                                  ScalarEvolution  **pointer_scalar_evolution_in,
                                  TensorEvolution  **pointer_tensor_evolution_in,
                                  int caseIndex_in)
:NDMaterial(tag, ND_TAG_NewTemplate3Dep), caseIndex(caseIndex_in)
{ 
    if ( pointer_material_parameter_in )
      pointer_material_parameter = pointer_material_parameter_in->newObj();
    else {
      opserr << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the input parameters. " << endln;
      exit(1);
    }

    if ( pointer_elastic_state_in )
      pointer_elastic_state = pointer_elastic_state_in->newObj();     
    else{
      opserr << "NewTemplate3Dep:: NewTemplate3Dep failed to get copy of elastic material. " << endln;
      exit(1);
    }

    if ( pointer_yield_function_in )
      pointer_yield_function = pointer_yield_function_in->newObj();
    else {
      opserr << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the yield function. " << endln;
      exit(1);
    }

    if ( pointer_plastic_flow_in )
      pointer_plastic_flow = pointer_plastic_flow_in->newObj();
    else {
      opserr << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the plastic flow. " << endln;
      exit(1);
    }

    // Scalar (isotropic) Evolution
    if ( pointer_material_parameter_in->getNum_Internal_Scalar() > 0 ) {
      pointer_scalar_evolution = new ScalarEvolution* [pointer_material_parameter_in->getNum_Internal_Scalar()];
      for (int i = 0; i < pointer_material_parameter_in->getNum_Internal_Scalar(); i++)
        pointer_scalar_evolution[i] = (pointer_scalar_evolution_in[i])->newObj();
    }
    else
      pointer_scalar_evolution = NULL;
    
    // Tensor (kinematic) Evolution
    if ( pointer_material_parameter_in->getNum_Internal_Tensor() > 0 ) {
      pointer_tensor_evolution = new TensorEvolution* [pointer_material_parameter_in->getNum_Internal_Tensor()];
      for (int i = 0; i < pointer_material_parameter_in->getNum_Internal_Tensor(); i++)
        pointer_tensor_evolution[i] = (pointer_tensor_evolution_in[i])->newObj();
    }
    else
      pointer_tensor_evolution = NULL;

    int err;
    err = this->revertToStart();
}

// Constructor
//================================================================================
NewTemplate3Dep::NewTemplate3Dep( int tag,
                                  MaterialParameter *pointer_material_parameter_in,
                                  ElasticState      *pointer_elastic_state_in,
                                  YieldFunction     *pointer_yield_function_in ,        
                                  PlasticFlow       *pointer_plastic_flow_in,
                                  TensorEvolution **pointer_tensor_evolution_in,
                                  int caseIndex_in)
:NDMaterial(tag, ND_TAG_NewTemplate3Dep), caseIndex(caseIndex_in)
{ 
    if ( pointer_material_parameter_in )
      pointer_material_parameter = pointer_material_parameter_in->newObj();
    else {
      opserr << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the input parameters. " << endln;
      exit(1);
    }

    if ( pointer_elastic_state_in )
      pointer_elastic_state = pointer_elastic_state_in->newObj();     
    else{
      opserr << "NewTemplate3Dep:: NewTemplate3Dep failed to get copy of elastic material. " << endln;
      exit(1);
    }

    if ( pointer_yield_function_in )
      pointer_yield_function = pointer_yield_function_in->newObj();
    else {
      opserr << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the yield function. " << endln;
      exit(1);
    }

    if ( pointer_plastic_flow_in )
       pointer_plastic_flow = pointer_plastic_flow_in->newObj();
    else {
      opserr << "NewTemplate3Dep:: NewTemplate3Dep failed to construct the plastic flow. " << endln;
      exit(1);
    }

    // Scalar (isotropic) Evolution
    pointer_scalar_evolution = NULL;
    
    // Tensor (kenimatic) Evolution
    if ( pointer_material_parameter_in->getNum_Internal_Tensor() > 0 ) {
      pointer_tensor_evolution = new TensorEvolution* [pointer_material_parameter_in->getNum_Internal_Tensor()]; 
      for (int i = 0; i < pointer_material_parameter_in->getNum_Internal_Tensor(); i++)
        pointer_tensor_evolution[i] = pointer_tensor_evolution_in[i]->newObj();
    }
    else
      pointer_tensor_evolution = NULL;

    int err;
    err = this->revertToStart();
}

// Constructor
//================================================================================
NewTemplate3Dep::NewTemplate3Dep()
: NDMaterial(0, ND_TAG_NewTemplate3Dep),
  pointer_material_parameter(NULL), 
  pointer_elastic_state(NULL), 
  pointer_yield_function(NULL), 
  pointer_plastic_flow(NULL), 
  pointer_scalar_evolution(NULL), 
  pointer_tensor_evolution(NULL), 
  caseIndex(0)
{
    int err; 
    err = this->revertToStart();
}

// Destructor
//================================================================================
NewTemplate3Dep::~NewTemplate3Dep()
{
    if (pointer_elastic_state)
      delete pointer_elastic_state;

    for (int i = 0; i < pointer_material_parameter->getNum_Internal_Scalar(); i++) {
      if (pointer_scalar_evolution[i])
        delete pointer_scalar_evolution[i];
    }
    if (pointer_scalar_evolution)
      delete [] pointer_scalar_evolution;         

    for (int j = 0; j < pointer_material_parameter->getNum_Internal_Tensor(); j++) {
      if (pointer_tensor_evolution[j])
        delete pointer_tensor_evolution[j];
    }
    if (pointer_tensor_evolution)
       delete [] pointer_tensor_evolution;

    if (pointer_yield_function)
       delete pointer_yield_function;

    if (pointer_plastic_flow)
       delete pointer_plastic_flow;      
    
    if (pointer_material_parameter)
       delete pointer_material_parameter;
}

// For Matrix based elements
int NewTemplate3Dep::setTrialStrain (const Vector &v)
{
	straintensor temp;
    
    temp.val(1,1) = v(0);  
    temp.val(2,2) = v(1);
    temp.val(3,3) = v(2);
    temp.val(1,2) = 0.5 * v(3);  
    temp.val(2,1) = 0.5 * v(3);  
    temp.val(3,1) = 0.5 * v(4);  
    temp.val(1,3) = 0.5 * v(4);
    temp.val(2,3) = 0.5 * v(5);  
    temp.val(3,2) = 0.5 * v(5);
    
	return this->setTrialStrainIncr(temp - getStrainTensor());
}

// For Matrix based elements
int NewTemplate3Dep::setTrialStrain (const Vector &v, const Vector &r)
{
	return this->setTrialStrainIncr(v);;
}

// For Matrix based elements
int NewTemplate3Dep::setTrialStrainIncr (const Vector &v)
{
	straintensor temp;
    
    temp.val(1,1) = v(0);  
    temp.val(2,2) = v(1);
    temp.val(3,3) = v(2);
    temp.val(1,2) = 0.5 * v(3);  
    temp.val(2,1) = 0.5 * v(3);  
    temp.val(3,1) = 0.5 * v(4);  
    temp.val(1,3) = 0.5 * v(4);
    temp.val(2,3) = 0.5 * v(5);  
    temp.val(3,2) = 0.5 * v(5);
	
    return this->setTrialStrainIncr(temp);
}

// For Matrix based elements
int NewTemplate3Dep::setTrialStrainIncr (const Vector &v, const Vector &r)
{
	return this->setTrialStrainIncr(v);
}

// For Matrix based elements
const Matrix& NewTemplate3Dep::getTangent (void)
{
   D(0,0) = Stiffness.cval(1,1,1,1);
   D(0,1) = Stiffness.cval(1,1,2,2);
   D(0,2) = Stiffness.cval(1,1,3,3);      
   D(0,3) = Stiffness.cval(1,1,1,2);
   D(0,4) = Stiffness.cval(1,1,1,3);
   D(0,5) = Stiffness.cval(1,1,2,3);      
    
   D(1,0) = Stiffness.cval(2,2,1,1);
   D(1,1) = Stiffness.cval(2,2,2,2);
   D(1,2) = Stiffness.cval(2,2,3,3);      
   D(1,3) = Stiffness.cval(2,2,1,2);
   D(1,4) = Stiffness.cval(2,2,1,3);
   D(1,5) = Stiffness.cval(2,2,2,3);            
    
   D(2,0) = Stiffness.cval(3,3,1,1);
   D(2,1) = Stiffness.cval(3,3,2,2);
   D(2,2) = Stiffness.cval(3,3,3,3);      
   D(2,3) = Stiffness.cval(3,3,1,2);
   D(2,4) = Stiffness.cval(3,3,1,3);
   D(2,5) = Stiffness.cval(3,3,2,3);                  
    
   D(3,0) = Stiffness.cval(1,2,1,1);
   D(3,1) = Stiffness.cval(1,2,2,2);
   D(3,2) = Stiffness.cval(1,2,3,3);      
   D(3,3) = Stiffness.cval(1,2,1,2);
   D(3,4) = Stiffness.cval(1,2,1,3);
   D(3,5) = Stiffness.cval(1,2,2,3);                        
    
   D(4,0) = Stiffness.cval(1,3,1,1);
   D(4,1) = Stiffness.cval(1,3,2,2);
   D(4,2) = Stiffness.cval(1,3,3,3);      
   D(4,3) = Stiffness.cval(1,3,1,2);
   D(4,4) = Stiffness.cval(1,3,1,3);
   D(4,5) = Stiffness.cval(1,3,2,3);                              
    
   D(5,0) = Stiffness.cval(2,3,1,1);
   D(5,1) = Stiffness.cval(2,3,2,2);
   D(5,2) = Stiffness.cval(2,3,3,3);      
   D(5,3) = Stiffness.cval(2,3,1,2);
   D(5,4) = Stiffness.cval(2,3,1,3);
   D(5,5) = Stiffness.cval(2,3,2,3);    

   return D;
}

// For Matrix based elements
const Vector& NewTemplate3Dep::getStress (void)
{
   sigma(0) = TrialStress.cval(1,1);
   sigma(1) = TrialStress.cval(2,2);
   sigma(2) = TrialStress.cval(3,3);
   sigma(3) = TrialStress.cval(1,2);
   sigma(4) = TrialStress.cval(1,3);
   sigma(5) = TrialStress.cval(2,3);

   return sigma;
}

// For Matrix based elements
const Vector& NewTemplate3Dep::getStrain (void)
{
   epsilon(0) = TrialStrain.cval(1,1);
   epsilon(1) = TrialStrain.cval(2,2);
   epsilon(2) = TrialStrain.cval(3,3);
   epsilon(3) = TrialStrain.cval(1,2) + TrialStrain.cval(2,1);
   epsilon(4) = TrialStrain.cval(1,3) + TrialStrain.cval(3,1);
   epsilon(5) = TrialStrain.cval(2,3) + TrialStrain.cval(3,2);    
    
   return epsilon;
}

//================================================================================
int NewTemplate3Dep::setTrialStrain(const Tensor& v)
{       
    return setTrialStrainIncr( v - getStrainTensor() ); 
}


//================================================================================
int NewTemplate3Dep::setTrialStrain(const Tensor& v, const Tensor& r)
{
    return setTrialStrainIncr( v - getStrainTensor() );
}

//================================================================================
int NewTemplate3Dep::setTrialStrainIncr(const Tensor& v)
{
   TrialStrain = v + getStrainTensor();
   
   switch(caseIndex) {
     
     case (0):
       return ForwardEuler(v);

     case (1):
       return SemiImplicit(v);
       
     case (2):
       return BackwardEuler(v);
       
     default:
       return 1;
    }
}

//================================================================================
int NewTemplate3Dep::setTrialStrainIncr(const Tensor& v, const Tensor& r)
{
    return setTrialStrainIncr(v);
}

//================================================================================
double NewTemplate3Dep::getRho(void)
{   
    double rho = 0.0;
    if (pointer_material_parameter->getNum_Material_Parameter() > 0)
        rho = pointer_material_parameter->getMaterial_Parameter(0);
    else {
      opserr << "Error!! NewTemplate3Dep:: number of input parameter for material constants less than 1. " << endln;
      opserr << "Remind: NewTemplate3Dep:: the 1st material constant is the density. " << endln;
      exit(1);
    }
    
    return rho;
}

//================================================================================
const BJtensor& NewTemplate3Dep::getTangentTensor(void)
{
    return Stiffness;
}

//================================================================================
const stresstensor&  NewTemplate3Dep::getStressTensor(void)
{
    return TrialStress;
}


//================================================================================
const straintensor& NewTemplate3Dep::getStrainTensor(void)
{
    return TrialStrain;
}

//================================================================================
const straintensor& NewTemplate3Dep::getPlasticStrainTensor(void)
{
    return TrialPlastic_Strain;
}


//================================================================================
int NewTemplate3Dep::commitState(void)
{
    int err = 0;
        
    //err += pointer_elastic_state->commitState();
    
    CommitStress.Initialize(TrialStress);
    CommitStrain.Initialize(TrialStrain);
    
    CommitPlastic_Strain.Initialize(TrialPlastic_Strain);
    
    return err;
}

//================================================================================
int NewTemplate3Dep::revertToLastCommit(void)
{
    int err = 0;
    
    TrialStress.Initialize(CommitStress);
    TrialStrain.Initialize(CommitStrain);
    
    TrialPlastic_Strain.Initialize(CommitPlastic_Strain);
    
    return err;
}

//================================================================================
int NewTemplate3Dep::revertToStart(void)
{
    int err = 0;
    
    CommitStress = pointer_elastic_state->getStress();
    CommitStrain = pointer_elastic_state->getStrain();
    
    CommitPlastic_Strain.Initialize(ZeroStrain);

    TrialStress.Initialize(CommitStress);
    TrialStrain.Initialize(CommitStrain);
    
    TrialPlastic_Strain.Initialize(ZeroStrain);

    Stiffness = pointer_elastic_state->getElasticStiffness(*pointer_material_parameter);
    
    return err;
}

//================================================================================
NDMaterial * NewTemplate3Dep::getCopy(void)
{
   NDMaterial* tmp = new NewTemplate3Dep(this->getTag(),
                                         this->pointer_material_parameter,
                                         this->pointer_elastic_state,
                                         this->pointer_yield_function,
                                         this->pointer_plastic_flow,
                                         this->pointer_scalar_evolution,
                                         this->pointer_tensor_evolution,
                                         this->caseIndex );
    return tmp;
}


//================================================================================
NDMaterial * NewTemplate3Dep::getCopy(const char *code)
{
    if (strcmp(code,"ThreeDimensional") == 0) {
       NewTemplate3Dep* tmp = new NewTemplate3Dep( this->getTag(),
                                                   this->pointer_material_parameter,
                                                   this->pointer_elastic_state,
                                                   this->pointer_yield_function,
                                                   this->pointer_plastic_flow,
                                                   this->pointer_scalar_evolution,
                                                   this->pointer_tensor_evolution,
                                                   this->caseIndex );
       return tmp;
    }
    else {
      opserr << "NewTemplate3Dep::getCopy failed to get model: " <<  code << endln;
      exit(1);
    }
    return 0;
}

//================================================================================
const char *NewTemplate3Dep::getType(void) const
{
    return "ThreeDimensional";
}

//================================================================================
int NewTemplate3Dep::sendSelf(int commitTag, Channel &theChannel)
{
    // Not yet implemented
    return 0;
}

//================================================================================
int NewTemplate3Dep::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // Not yet implemented
    return 0;
}

//================================================================================
void NewTemplate3Dep::Print(OPS_Stream& s, int flag)
{
     s << (*this);
}

//================================================================================
int NewTemplate3Dep::ForwardEuler(const straintensor& strain_incr)
{
    straintensor start_strain;
    stresstensor start_stress;
    stresstensor stress_incr;
    stresstensor Intersection_strain;
    stresstensor Intersection_stress;
    stresstensor elastic_predictor_stress;
    BJtensor Ee;
    straintensor  incr_strain;
    int err = 0;

    double f_start = 0.0;
    double f_pred  = 0.0;
    
    double intersection_factor = 0.0;
    
    Ee = pointer_elastic_state->getElasticStiffness(*pointer_material_parameter);
    
    start_stress = getStressTensor();
    start_strain = getStrainTensor();
    
    Intersection_stress.Initialize(start_stress);
    
    // I had to use the one line incr_strain = strain_incr; Problem of BJTensor;
    incr_strain.Initialize(strain_incr);
    stress_incr = Ee("ijpq") * incr_strain("pq");
    stress_incr.null_indices();

    elastic_predictor_stress = start_stress + stress_incr;
   
    f_start = pointer_yield_function->YieldFunctionValue( start_stress, *pointer_material_parameter );
    f_pred =  pointer_yield_function->YieldFunctionValue( elastic_predictor_stress, *pointer_material_parameter );
    
    // If Elastic
    if ( (f_start <= 0.0 && f_pred <= FTOL) || f_start > f_pred ) {
        //TrialStrain = start_strain + strain_incr;
        TrialStress.Initialize(elastic_predictor_stress);

        Stiffness = Ee;
               
        // update elastic part
        err += pointer_elastic_state->setStress(TrialStress);
        err += pointer_elastic_state->setStrain(TrialStrain);       
                
        return err;
    }
    
    // If Elastic and then Elastic-Plastic
    if ( f_start < 0.0 )  {

        intersection_factor = zbrentstress( start_stress, elastic_predictor_stress, 0.0, 1.0, TOL );      
        
        Intersection_stress = yield_surface_cross( start_stress, elastic_predictor_stress, intersection_factor );
        Intersection_strain = start_strain + (incr_strain * intersection_factor);         
        
        stress_incr = elastic_predictor_stress - Intersection_stress;
        
        // Update elastic part
        err += pointer_elastic_state->setStress(Intersection_stress);      
        err += pointer_elastic_state->setStrain(Intersection_strain);
        
        Ee = pointer_elastic_state->getElasticStiffness(*pointer_material_parameter);
    }
    
    // If E-P Response,
    {
        double lower = 0.0;
        double Delta_lambda = 0.0;
        double hardMod  = 0.0;
        double h_s = 0.0;
        double xi_s = 0.0;
        stresstensor dFods;
        straintensor dQods;
        BJtensor Hq;
        BJtensor Hf;
        BJtensor Ep;
   
        stresstensor h_t;
        stresstensor xi_t;
        
        straintensor plastic_strain_incr;
        stresstensor ep_stress;
        
        // For better numerical performance
        //Intersection_stress = Intersection_stress *(1.0 - TOL);
                   
        dFods = pointer_yield_function->StressDerivative( Intersection_stress, *pointer_material_parameter );
        dQods = pointer_plastic_flow->PlasticFlowTensor( Intersection_stress, Intersection_strain, *pointer_material_parameter );

        // E_ijkl * R_kl
        Hq = Ee("ijkl") * dQods("kl");
        Hq.null_indices();
        
        // L_ij * E_ijkl
        Hf = dFods("ij") * Ee("ijkl");
        Hf.null_indices();
        
        // L_ij * E_ijkl * R_kl
        lower = ( Hf("ij") * dQods("ij") ).trace(); 
        int i;  
        // Evolution of scalar (isotropic) internal variables in yield function
        double Num_internal_scalar_in_yield_function = pointer_yield_function->getNumInternalScalar();
        for (i = 0; i < Num_internal_scalar_in_yield_function; i++) {
          h_s = pointer_scalar_evolution[i]->H( dQods, Intersection_stress, Intersection_strain, *pointer_material_parameter);
          xi_s = pointer_yield_function->InScalarDerivative( Intersection_stress, *pointer_material_parameter, i+1);
          hardMod += h_s * xi_s;
        }
        
        // Evolution of tensor (kinematic) internal variables in yield function
        double Num_internal_tensor_in_yield_function = pointer_yield_function->getNumInternalTensor();
        for (i = 0;  i < Num_internal_tensor_in_yield_function; i++) {
          h_t = pointer_tensor_evolution[i]->Hij( dQods, Intersection_stress, Intersection_strain, *pointer_material_parameter);
          xi_t = pointer_yield_function->InTensorDerivative( Intersection_stress, *pointer_material_parameter, i+1);
          hardMod += ( h_t("mn") * xi_t("mn") ).trace();
        }

        lower -= hardMod;        
                
        // L_ij * E_ijkl * d e_kl ( true ep strain increment)
        Delta_lambda = ( dFods("ij") * stress_incr("ij") ).trace();
        
        if (lower != 0.0)
          Delta_lambda /= lower;              
         
        if (Delta_lambda < 0.0)  
          Delta_lambda = 0.0;

        // Plastic strain increment
        plastic_strain_incr = dQods * Delta_lambda;
        ep_stress = elastic_predictor_stress - (Hq * Delta_lambda);

        TrialPlastic_Strain = this->getPlasticStrainTensor() + plastic_strain_incr;
        TrialStress = ep_stress;
        
        // To obtain Eep
        Ep = Hq("pq") * Hf("mn");  
        Ep.null_indices();        	
        
        Ep = Ep * (1.0/lower);
        if ( Delta_lambda > 0.0 )
        	Stiffness = Ee - Ep;
        else
        	Stiffness = Ee;

        // Update internal scalar variables
        double dS= 0.0;
        double S = 0.0;
        int Num_internal_scalar = pointer_material_parameter->getNum_Internal_Scalar();
        for (i = 0; i < Num_internal_scalar; i++) {
          dS = ( pointer_scalar_evolution[i]->H(dQods, Intersection_stress, Intersection_strain, *pointer_material_parameter) ) *Delta_lambda;
          S = pointer_material_parameter->getInternal_Scalar(i);
          err += pointer_material_parameter->setInternal_Scalar(i, S + dS );
        }
        
        // Update internal tensor variables
        stresstensor dT;
        stresstensor T;
        int Num_internal_tensor = pointer_material_parameter->getNum_Internal_Tensor();
        for (i = 0; i < Num_internal_tensor; i++) {
          dT = pointer_tensor_evolution[i]->Hij(dQods, Intersection_stress, Intersection_strain, *pointer_material_parameter) *Delta_lambda; 
          T = pointer_material_parameter->getInternal_Tensor(i);
          err += pointer_material_parameter->setInternal_Tensor(i, T + dT );
        }
        
        // Update elastic part
        err += pointer_elastic_state->setStrain(TrialStrain);
        err += pointer_elastic_state->setStress(TrialStress);       
    }

    return err;
}

//================================================================================
int NewTemplate3Dep::SemiImplicit(const straintensor& strain_incr)
{
    double YieldFun = 0.0;

    straintensor start_strain;
    stresstensor start_stress;
    stresstensor stress_incr;
    BJtensor Ee;
 
    int err = 0;

    Ee = pointer_elastic_state->getElasticStiffness(*pointer_material_parameter);
    
    start_stress = getStressTensor();
    start_strain = getStrainTensor();
    
    // I had to use the one line incr_strain = strain_incr; Problem of BJTensor;
    straintensor  incr_strain = strain_incr;
    stress_incr = Ee("ijpq") * incr_strain("pq");    
    stress_incr.null_indices();

    TrialPlastic_Strain = this->getPlasticStrainTensor();
    TrialStress = start_stress + stress_incr;

    YieldFun = pointer_yield_function->YieldFunctionValue(TrialStress, *pointer_material_parameter);
    
    //opserr << "YieldFun = " << YieldFun << endln;
    
    if ( YieldFun <= FTOL ) {      // If Elastic
        err += pointer_elastic_state->setStress(TrialStress);
        err += pointer_elastic_state->setStrain(TrialStrain);       
        Stiffness = Ee;
    }
    else {      // If Elastic-Plastic
        double Delta_lambda  = 0.0;
        double d2_lambda  = 0.0;
        int	 iter_counter = 0;
        double lower = 0.0;
        double hardMod  = 0.0;
        double h_s = 0.0;
        double xi_s = 0.0;
        stresstensor dFods;
        straintensor dQods;
        BJtensor Hf;
        BJtensor Hq;
        BJtensor Ep;
   
        stresstensor h_t;
        stresstensor xi_t;
               
        dQods = pointer_plastic_flow->PlasticFlowTensor( start_stress, start_strain, *pointer_material_parameter );

        Hq = Ee("ijkl") * dQods("kl");
          Hq.null_indices();
        int i;
        // Evolution of scalar (isotropic) internal variables in yield function
        double Num_internal_scalar_in_yield_function = pointer_yield_function->getNumInternalScalar();
        for (i = 0; i < Num_internal_scalar_in_yield_function; i++) {
          h_s = pointer_scalar_evolution[i]->H( dQods, start_stress, start_strain, *pointer_material_parameter);
          xi_s = pointer_yield_function->InScalarDerivative( TrialStress, *pointer_material_parameter, i+1);
          hardMod += h_s * xi_s;
        }        
        // Evolution of tensor (kinematic) internal variables in yield function
        double Num_internal_tensor_in_yield_function = pointer_yield_function->getNumInternalTensor();
        for (i = 0;  i < Num_internal_tensor_in_yield_function; i++) {
          h_t = pointer_tensor_evolution[i]->Hij( dQods, start_stress, start_strain, *pointer_material_parameter);
          xi_t = pointer_yield_function->InTensorDerivative( TrialStress, *pointer_material_parameter, i+1);
          hardMod += ( h_t("mn") * xi_t("mn") ).trace();
        }

        // ################## Beginning of do-while ########################
        do {           
          dFods = pointer_yield_function->StressDerivative( TrialStress, *pointer_material_parameter );
        
          Hf = dFods("ij") * Ee("ijkl");
          Hf.null_indices();        

          lower = ( Hf("ij") * dQods("ij") ).trace();          
          lower = lower - hardMod;        
                
          d2_lambda = YieldFun / lower;              
        
          // Update stress
          TrialStress -= (Hq *d2_lambda);

          // Update internal scalar variables
          double dS= 0.0;
          double S = 0.0;
		  int i;
          int Num_internal_scalar = pointer_material_parameter->getNum_Internal_Scalar();
          for (i = 0; i < Num_internal_scalar; i++) {
            dS = ( pointer_scalar_evolution[i]->H(dQods, start_stress, start_strain, *pointer_material_parameter) ) *d2_lambda;
            S = pointer_material_parameter->getInternal_Scalar(i);
            err += pointer_material_parameter->setInternal_Scalar(i, S + dS );
          }        
          // Update internal tensor variables
          stresstensor dT;
          stresstensor T;
          int Num_internal_tensor = pointer_material_parameter->getNum_Internal_Tensor();
          for (i = 0; i < Num_internal_tensor; i++) {
            dT = pointer_tensor_evolution[i]->Hij(dQods, start_stress, start_strain, *pointer_material_parameter) *d2_lambda; 
            T = pointer_material_parameter->getInternal_Tensor(i);
            err += pointer_material_parameter->setInternal_Tensor(i, T + dT );
          }

          // Update elastic part
          err += pointer_elastic_state->setStress(TrialStress);

          // Update Delta_lambda
          Delta_lambda += d2_lambda;

          // Update iter_counter
          iter_counter++;
          
          // Update Yield Function
          YieldFun = pointer_yield_function->YieldFunctionValue(TrialStress, *pointer_material_parameter);
          //opserr << "F = " << YieldFun << endln;

          //if (iter_counter == ITMAX)
          //  opserr << "Warning! The iteration number in the semi-implicit algorithm reaches to " << ITMAX << endln;

        } while (YieldFun > FTOL && iter_counter < ITMAX);
        // ################## End of do-while ########################
        
        if (Delta_lambda < 0.0)
          Delta_lambda = 0.0;
       
        // Return algorithmic stiffness tensor
        Ep = Hq("pq") * Hf("mn");  
          Ep.null_indices();   
        Ep = Ep * (1.0/lower);
        if ( Delta_lambda > 0.0 )
        	Stiffness = Ee - Ep;
        else
        	Stiffness = Ee;

        TrialPlastic_Strain += (dQods *Delta_lambda);
        
        err += pointer_elastic_state->setStrain(TrialStrain);
    }

    return err;
}

//================================================================================
int NewTemplate3Dep::BackwardEuler(const straintensor& strain_incr)
{
	// Need Work Here!
    opserr << "BackwardEuler is not yet implemented!" << endln;
	return 0;
}


// Trying to find intersection point according to M. Crisfield's book
// "Non-linear Finite Element Analysis of Solids and Structures "  Chp 6.6.1 pp168.
//================================================================================
stresstensor NewTemplate3Dep::yield_surface_cross(const stresstensor & start_stress,
                                                  const stresstensor & end_stress, double a)
{
    //double a = zbrentstress( start_stress, end_stress, 0.0, 1.0, TOL );

    stresstensor delta_stress = end_stress - start_stress;
    stresstensor intersection_stress = start_stress + delta_stress * a;

    return intersection_stress;
}

// Routine used by yield_surface_cross to find the stresstensor at cross point
//================================================================================
double NewTemplate3Dep::zbrentstress(const stresstensor& start_stress,
                                  const stresstensor& end_stress,
                                  double x1, double x2, double tol) const
{
  double EPS = d_macheps();

  int iter;
  double a = x1;
  double b = x2;
  double c = 0.0;
  double d = 0.0;
  double e = 0.0;
  double min1 = 0.0;
  double min2 = 0.0;
  double fc = 0.0;
  double p = 0.0;
  double q = 0.0;
  double r = 0.0;
  double s = 0.0;
  double tol1 = 0.0;
  double xm = 0.0;  
   
  double fa = func(start_stress, end_stress, *pointer_material_parameter, a);
  double fb = func(start_stress, end_stress, *pointer_material_parameter, b);
 
  if ( (fb * fa) > 0.0) {
      opserr << "\a\n Root must be bracketed in ZBRENTstress " << endln;
      exit(1);
  }
  
  fc = fb;
  for ( iter = 1; iter <= ISMAX; iter++ ) {
      if ( (fb * fc) > 0.0) {
          c = a;
          fc = fa;
          e = d = b - a;
      }
      if ( fabs(fc) < fabs(fb) ) { 
          a = b;
          b = c;
          c = a;
          fa = fb;
          fb = fc;
          fc = fa;
      }
      tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
      xm = 0.5 * (c - b);
      if ( fabs(xm) <= tol1 || fb == 0.0 ) 
        return b;
        
      if ( fabs(e) >= tol1 && fabs(fa) > fabs(fb) ) {
          s = fb / fa;
          if (a == c) {
              p = 2.0 * xm * s;
              q = 1.0 - s;
          }
          else {
              q = fa / fc;
              r = fb / fc;
              p = s * ( 2.0 * xm * q * (q - r) - (b - a) * (r - 1.0) );
              q = (q - 1.0) * (r - 1.0) * (s - 1.0);
          }
          if (p > 0.0)  
            q = -q;
          p = fabs(p);
          min1 = 3.0 * xm * q - fabs(tol1*q);
          min2 = fabs(e*q);
          if (2.0*p < (min1 < min2 ? min1 : min2)) {
              e = d;
              d = p/q;
          }
          else {
              d = xm;
              e = d;
          }
        }
        else {
          d = xm;
          e = d;
        }
      a = b;
      fa = fb;
      if (fabs(d) > tol1)
        b += d;
      else
        b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
      fb = func(start_stress, end_stress, *pointer_material_parameter, b);
  }
  
  return 0.0;
}

//================================================================================
double NewTemplate3Dep::func(const stresstensor& start_stress,
                          const stresstensor& end_stress,
                          const MaterialParameter& pointer_material_parameter, 
                          double alfa ) const
{     
    stresstensor alfa_stress = ( start_stress * (1.0 - alfa) ) + ( end_stress * alfa );

    double f = pointer_yield_function->YieldFunctionValue( alfa_stress, pointer_material_parameter );
   
    return f;
}

#endif

