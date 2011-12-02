/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elaso-plastic constitutive model       #
#                    implementation                                              #
# CLASS:             DPEPState (the base class for all Elasto-plastic state)     #
#                                                                                #
# VERSION:                                                                       #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )    #
# TARGET OS:         DOS || UNIX || . . .                                        #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                 #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                 #
#                                                                                #
#                                                                                #
# DATE:              08-03-2000                                                  #
# UPDATE HISTORY:                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
# SHORT EXPLANATION: This class is used to hold all state parameters and internal#
#                    variables in an elasto-plastic constitutive model!          #
#                                                                                #
//================================================================================
*/

#ifndef EPState_CPP
#define EPState_CPP

#include "EPState.h"
#include <G3Globals.h>

//================================================================================
//Normal Constructor 1
//================================================================================

EPState::EPState(double               Eod,
                 double               Ed,
                 double               nu,
                 double               rho,
                 const stresstensor  &stressp,       
                 const straintensor  &strainp, 
                 const straintensor  &Estrainp,
                 const straintensor  &Pstrainp,
                 const straintensor  &dEstrainp,
                 const straintensor  &dPstrainp,
	         int                  NScalarp,
		 const double       * Scalarp,
	         int                  NTensorp,
	         const stresstensor * Tensorp,
	         const tensor       & Eepp, 
    	         const stresstensor & Stress_commitp,
    	         const straintensor & Strain_commitp,	  
    	         const double       * Scalar_commitp,
    	         const stresstensor * Tensor_commitp, 
    	         const tensor       & Eep_commitp,
    	         const stresstensor & Stress_initp,    
    	         const straintensor & Strain_initp,	 
    	         const double       * Scalar_initp,
    	         const stresstensor * Tensor_initp,
    	         const tensor       & Eep_initp, 
                 bool               Convergedp) 
: Eo(Eod), E_Young(Ed), nu_Poisson(nu), rho_mass_density(rho), CurrentStress(stressp),
  CurrentStrain(strainp), ElasticStrain(Estrainp), PlasticStrain(Pstrainp), 
  dElasticStrain(dEstrainp), dPlasticStrain(dPstrainp), Eep(Eepp), 
  Stress_commit(Stress_commitp), Strain_commit(Strain_commitp),
  Eep_commit(Eep_commitp), Stress_init(Stress_initp), Strain_init(Strain_initp), 
  Eep_init(Eep_initp), Converged (Convergedp)  
{

      Eo               = Eod;	        
      E_Young          = Ed;	        
      nu_Poisson       = nu;	     
      rho_mass_density = rho; 
      
      NScalarVar = NScalarp;
      //ScalarVar = new double[ NScalarVar ]; 
      //if ( !ScalarVar ) {
      if ( ( !Scalarp ) && (NScalarVar != 0) ) {
         g3ErrorHandler->warning("EPState::EPState   No initial values for scalar hardening vars, set to aero");
         //::exit(1);  
      }
      else {
         for (int i = 0; i < NScalarVar; i++) {
            //cout << Scalarp[i] << endln; 
        	 ScalarVar[i] = Scalarp[i];
        	 ScalarVar_commit[i] = Scalar_commitp[i];
        	 ScalarVar_init[i] = Scalar_initp[i];
         }
      }

      NTensorVar = NTensorp;
      //TensorVar = new stresstensor[ NTensorVar ];
      //if ( !TensorVar ) {
      //   g3ErrorHandler->fatal("EPState::EPState insufficient memory for Tensor hardening vars");
      //   ::exit(1);  
      //}

      if ( (!Tensorp ) && ( NTensorVar)) {
         //ScalarVar = new double[ NScalarVar ]; 
         //if ( !ScalarVar ) {
         g3ErrorHandler->warning("EPState::EPState   No initial values for tensorial hardening vars, set to zero");
         //::exit(1);  
      }
      else {      
         for (int i = 0; i < NTensorVar; i++) {
         	 //cout << Tensorp[i];
         	 TensorVar[i] = Tensorp[i];
         	 TensorVar_commit[i] = Tensor_commitp[i];
         	 TensorVar_init[i] = Tensor_initp[i];
         	 //cout << TensorVar[i];
         	 //TensorVar[i].null_indices();
         }
      }
}

//================================================================================
//Normal Constructor 11
//================================================================================

EPState::EPState(double             Eod,
                 double             Ed,
                 double             nu,
                 double             rho,
                 const stresstensor &stressp,       
                 const straintensor &strainp, 
                 const straintensor &Estrainp,
                 const straintensor &Pstrainp,
	         int                NScalarp,
		 const double     * Scalarp,
	         int                NTensorp,
	         const tensor     * Tensorp )	  
: CurrentStress(stressp), CurrentStrain(strainp), ElasticStrain(Estrainp), 
  PlasticStrain(Pstrainp), Stress_commit(stressp), Strain_commit(strainp), 
  Stress_init(stressp), Strain_init(strainp)  
{
      Eo               = Eod;	        
      E_Young          = Ed;	        
      nu_Poisson       = nu;	     
      rho_mass_density = rho; 
      
      CurrentStress    = stressp;
      //cout << "stressp " << stressp;
      //CurrentStress.null_indices();
      CurrentStrain    =  strainp;
      ElasticStrain    =  Estrainp;
      PlasticStrain    =  Pstrainp;
      //Eep = Eepp;
      //cout << "strainp " << strainp;
      //CurrentStrain.null_indices();

      NScalarVar = NScalarp;
      //ScalarVar = new double[ NScalarVar ]; 
      //if ( !ScalarVar ) {
      if ( ( !Scalarp ) && (NScalarVar != 0) ) {
         g3ErrorHandler->warning("EPState::EPState   No initial values for scalar hardening vars, set to zero");
         //::exit(1);  
      }
      else {
         for (int i = 0; i < NScalarVar; i++) {
            //cout << Scalarp[i] << endln; 
       	 ScalarVar[i] = Scalarp[i];
       	 ScalarVar_commit[i] = Scalarp[i];
       	 ScalarVar_init[i] = Scalarp[i];
         }
      }

      NTensorVar = NTensorp;
      //TensorVar = new stresstensor[ NTensorVar ];
      //if ( !TensorVar ) {
      //   g3ErrorHandler->fatal("EPState::EPState insufficient memory for Tensor hardening vars");
      //   ::exit(1);  
      //}

      if ( ( !Tensorp ) && ( NTensorVar != 0 ) ) {
         //ScalarVar = new double[ NScalarVar ]; 
         //if ( !ScalarVar ) {
         g3ErrorHandler->warning("EPState::EPState   No initial values for tensorial hardening vars, set to zero");
         //::exit(1);  
      }
      else {       
         for (int i = 0; i < NTensorVar; i++) {
       	 //cout << Tensorp[i];
       	 TensorVar[i] = Tensorp[i];
       	 TensorVar_commit[i] = Tensorp[i];
       	 TensorVar_init[i] = Tensorp[i];
       	 //cout << TensorVar[i];
       	 //TensorVar[i].null_indices();
         }
     }

     Converged = false;
}

//================================================================================
//Normal Constructor 2
//================================================================================

EPState::EPState(double             Eod,
                 double             Ed,
                 double             nu,
                 double             rho,
	         int                NScalarp,
	         const double     * Scalarp,
	         int                NTensorp,
	         const tensor     * Tensorp ) {

      Eo               = Eod;	        
      E_Young          = Ed;	        
      nu_Poisson       = nu;	     
      rho_mass_density = rho; 

      //CurrentStress    = stressp;
      //cout << "CurrentStress " << CurrentStress;
      //CurrentStress.null_indices();
      //CurrentStrain    =  strainp;
      //ElasticStrain    =  Estrainp;
      //PlasticStrain    =  Pstrainp;
      //dElasticStrain   =  dEstrainp;
      //dPlasticStrain   =  dPstrainp;
      
      Eep = tensor( 4, def_dim_4, 0.0 ); // need to be initialized as 4th order tensor
      //cout << "strainp " << strainp;
      //CurrentStrain.null_indices();

      NScalarVar  =  NScalarp;
      //ScalarVar = new double[ NScalarVar ]; 
      //if ( !ScalarVar ) {
      //   g3ErrorHandler->fatal("EPState::EPState insufficient memory for Scalar hardening vars");
      //   ::exit(1);  
      //}
      if ( (!Scalarp) && (NScalarVar != 0) ) {
         g3ErrorHandler->warning("EPState::EPState   No initial values for scalar hardening vars, set to zero");
         //::exit(1);  
      }
      else {
         for (int i = 0; i < NScalarVar; i++) {
            //cout << Scalarp[i] << endln; 
        	 ScalarVar[i] = Scalarp[i];
        	 ScalarVar_commit[i] = Scalarp[i];
        	 ScalarVar_init[i] = Scalarp[i];
         }
      }

      NTensorVar = NTensorp;
      //TensorVar = new stresstensor[ NTensorVar ];
      //if ( !TensorVar ) {
      //   g3ErrorHandler->fatal("EPState::EPState insufficient memory for Tensor hardening vars");
      //   ::exit(1);  
      //}
      if ( (!Scalarp) && ( NTensorVar != 0 )) {
         g3ErrorHandler->warning("EPState::EPState   No initial values for tensorial hardening vars, set to zero");
         //::exit(1);  
      }
      else {
         for (int i = 0; i < NTensorVar; i++) {
        	 //cout << Tensorp[i];
        	 TensorVar[i] = Tensorp[i];
        	 TensorVar_commit[i] = Tensorp[i];
        	 TensorVar_init[i] = Tensorp[i];
        	 //cout << TensorVar[i];
        	 //TensorVar[i].null_indices();
         }
      }
      
      Converged = false;
}


//================================================================================
//Normal Constructor --no parameters
//================================================================================

EPState::EPState( ) {

      Eo               = 80000.0;
      E_Young          = 80000.0;
      nu_Poisson       = 0.3;	     
      rho_mass_density = 0.0; 
      Eep = tensor( 4, def_dim_4, 0.0 );

      NScalarVar = MaxNScalarVar;
      for (int i =0; i < NScalarVar; i++) 
         ScalarVar[i] = 0.0; 

      NTensorVar = MaxNTensorVar;
      //for (int i =0; i < NTensorVar, i++) 
      //   TensorVar[i] = stresstensor(0.0);
 
      Converged = false;

}

//================================================================================
//create a clone of itself
//================================================================================

EPState* EPState::newObj() {

      EPState * eps = new  EPState(this->getEo(),	       
      				   this->getE(),
      				   this->getnu(),
      				   this->getrho(),
                                   this->getStress(),
      				   this->getStrain(), 
      				   this->getElasticStrain(), 
      				   this->getPlasticStrain(), 
      				   this->getdElasticStrain(), 
      				   this->getdPlasticStrain(), 
      				   this->getNScalarVar(), 
      				   this->getScalarVar(),
      				   this->getNTensorVar(),
      				   this->getTensorVar(),
      				   this->getEep(), 
                                   this->getStress_commit(),
      				   this->getStrain_commit(), 
      				   this->getScalarVar_commit(),
      				   this->getTensorVar_commit(),
      				   this->getEep_commit(), 
                                   this->getStress_init(),
      				   this->getStrain_init(), 
      				   this->getScalarVar_init(),
      				   this->getTensorVar_init(),
      				   this->getEep_init(), 
      				   this->getConverged()
				   );
      return eps;
}      				 


//================================================================================
// Copy constructor
//================================================================================
EPState::EPState( const EPState &rhs ) {         

      Eo               = rhs.getEo();	        
      E_Young          = rhs.getE();	        
      nu_Poisson       = rhs.getnu();	     
      rho_mass_density = rhs.getrho(); 
      CurrentStress    = rhs.getStress();
      CurrentStrain    = rhs.getStrain();
      ElasticStrain    = rhs.getElasticStrain();
      PlasticStrain    = rhs.getPlasticStrain();
      dElasticStrain   = rhs.getdElasticStrain();
      dPlasticStrain   = rhs.getdPlasticStrain();
      
      Stress_commit = rhs.getStress_commit();
      Strain_commit = rhs.getStrain_commit();
      Stress_init   = rhs.getStress_init();
      Strain_init   = rhs.getStrain_init();
      //cout << Eep.rank() << " ";
      //Eep.printshort("before copy con ");
      Eep = rhs.getEep();
      Eep_commit = rhs.getEep_commit();
      Eep_init = rhs.getEep_init();
      //cout << Eep.rank() << " ";
      //Eep.printshort("after copy con ");

      NScalarVar  =  rhs.getNScalarVar();
      //ScalarVar = new double[ NScalarVar ]; 
      //if ( !ScalarVar ) {
      //   g3ErrorHandler->fatal("EPState::EPState insufficient memory for Scalar hardening vars");
      //   ::exit(1);  
      //}
      for (int i = 0; i < NScalarVar; i++) { 
	 ScalarVar[i] = rhs.ScalarVar[ i ];
	 ScalarVar_commit[i] = rhs.ScalarVar_commit[ i ];
	 ScalarVar_init[i] = rhs.ScalarVar_init[ i ];
      }
      NTensorVar = rhs.getNTensorVar();
      //TensorVar = new stresstensor[ NTensorVar ];
      //if ( !TensorVar ) {
      //   g3ErrorHandler->fatal("EPState::EPState insufficient memory for Tensor hardening vars");
      //   ::exit(1);  
      //}
      for (int i = 0; i < NTensorVar; i++) {
	 TensorVar[i] = rhs.TensorVar[ i ];
	 TensorVar_commit[i] = rhs.TensorVar_commit[ i ];
	 TensorVar_init[i] = rhs.TensorVar_init[ i ];
	 //cout << TensorVar[i];
 	 //TensorVar[i].null_indices();
      }
      
      Converged  =  rhs.getConverged();
     
}      				 


//================================================================================
//Overloading the assignment sign
//================================================================================
const EPState & EPState::operator=(const EPState &rhs ) { 
        
      if ( this != &rhs ) {
         Eo               = rhs.getEo();
         E_Young          = rhs.getE();	        
         nu_Poisson       = rhs.getnu();	     
         rho_mass_density = rhs.getrho(); 

         CurrentStress    = rhs.getStress();
         //cout << "Current stress " << CurrentStress;
         CurrentStrain    = rhs.getStrain();
         //cout << "strainp " << strainp;
         ElasticStrain    = rhs.getElasticStrain();
         PlasticStrain    = rhs.getPlasticStrain();
         dElasticStrain   = rhs.getdElasticStrain();
         dPlasticStrain   = rhs.getdPlasticStrain();

         Stress_commit = rhs.getStress_commit();
         Strain_commit = rhs.getStrain_commit();
         Stress_init   = rhs.getStress_init();
         Strain_init   = rhs.getStrain_init();
         
	 Eep              = rhs.getEep();
         Eep_commit = rhs.getEep_commit();
         Eep_init = rhs.getEep_init();
         
         NScalarVar  =  rhs.getNScalarVar();
         //ScalarVar = new double[ NScalarVar ]; 
         //if ( !ScalarVar ) {
         //   g3ErrorHandler->fatal("EPState::operator= insufficient memory for Scalar hardening vars");
         //   ::exit(1);  
         //}
         for (int i = 0; i < NScalarVar; i++) {
            ScalarVar[i] = rhs.ScalarVar[i];
            ScalarVar_commit[i] = rhs.ScalarVar_commit[i];
            ScalarVar_init[i] = rhs.ScalarVar_init[i];
         }
         
         NTensorVar = rhs.getNTensorVar();
         //TensorVar = new stresstensor[ NTensorVar ];
         //if ( !TensorVar ) {
         //   g3ErrorHandler->fatal("EPState::operator= insufficient memory for Tensor hardening vars");
         //   ::exit(1);  
         //}
         for (int i = 0; i < NTensorVar; i++) {
             TensorVar[i] = rhs.TensorVar[i];
             TensorVar_commit[i] = rhs.TensorVar_commit[i];
             TensorVar_init[i] = rhs.TensorVar_init[i];
             //TensorVar[i].null_indices();
         }

         Converged = rhs.getConverged();
      
      }			     
      
      return *this;
      
}      				 

//================================================================================
double EPState::getE() const {
      return E_Young; 
}

//================================================================================
double EPState::getEo() const {
      return Eo; 
}

//================================================================================
double EPState::getnu() const {
      return nu_Poisson; 
}

//================================================================================
double EPState::getrho() const {
      return rho_mass_density; 
};

//================================================================================
int EPState::getNScalarVar() const {
      return NScalarVar; 
}

//================================================================================
int EPState::getNTensorVar() const {
      return NTensorVar; 
}

//================================================================================
bool EPState::getConverged() const {
      return Converged; 
}
		      

//================================================================================
stresstensor EPState::getStress() const {

     return CurrentStress;

}

//================================================================================
stresstensor EPState::getStress_commit() const {

     return Stress_commit;

}

//================================================================================
stresstensor EPState::getStress_init() const {

     return Stress_init;

}

//================================================================================
stresstensor EPState::getIterativeStress() const {

     return IterativeStress;

}

//================================================================================
straintensor EPState::getStrain() const {

     return CurrentStrain;

}


//================================================================================
straintensor EPState::getStrain_commit() const {

     return Strain_commit;

}

//================================================================================
straintensor EPState::getStrain_init() const {

     return Strain_init;

}

//================================================================================
straintensor EPState::getElasticStrain() const {

     return ElasticStrain;

}

//================================================================================
straintensor EPState::getPlasticStrain() const {

     return PlasticStrain;

}
//================================================================================
straintensor EPState::getdElasticStrain() const {

     return dElasticStrain;

}


//================================================================================
straintensor EPState::getdPlasticStrain() const {

     return dPlasticStrain;

}


//================================================================================
tensor EPState::getEep() const {

     return Eep;

}

//================================================================================
tensor EPState::getEep_commit() const {

     return Eep_commit;

}

//================================================================================
tensor EPState::getEep_init() const {

     return Eep_init;

}


//================================================================================
void EPState::setE( double Ey ) { 
      E_Young = Ey; 
}


//================================================================================
void EPState::setStress(const stresstensor &newstress ) { 
      CurrentStress = newstress; 
}

//================================================================================
void EPState::setIterativeStress(const stresstensor &newstress ) { 
      IterativeStress = newstress; 
}


//================================================================================
void EPState::setStrain(const straintensor &newstrain ) {

      CurrentStrain = newstrain; 

}

//================================================================================
void EPState::setElasticStrain(const straintensor &newEstrain ) {

      ElasticStrain = newEstrain; 

}

//================================================================================
void EPState::setPlasticStrain(const straintensor &newPstrain ) {

      PlasticStrain = newPstrain; 

}


//================================================================================
void EPState::setdElasticStrain(const straintensor &newdEstrain ) {

      dElasticStrain = newdEstrain; 

}

//================================================================================
void EPState::setdPlasticStrain(const straintensor &newdPstrain ) {

      dPlasticStrain = newdPstrain; 

}

//================================================================================
void EPState::setEep(const tensor &newEep )  {

     Eep = newEep;

}

//================================================================================
void EPState::setConverged( bool b) {
     Converged = b;
}


//================================================================================
// Retrun the nth Scalar Var.... Starting from 1!!
//================================================================================
double EPState::getScalarVar( int WhichOne) const { 

      if ( WhichOne <= getNScalarVar() )  
         return ScalarVar[ WhichOne - 1 ]; 
      else 
      {
         g3ErrorHandler->fatal("EPState::getScalarVar Out of ScalarVar's range %d!", getNScalarVar() );
	 exit(1);
      }

}


//================================================================================
// Retrun the nth Tensor Var.... Starting from 1!!
//================================================================================
stresstensor EPState::getTensorVar(int WhichOne) const { 

      if ( WhichOne <= getNTensorVar() )  
         return TensorVar[ WhichOne - 1 ]; 
      else 
      {
         g3ErrorHandler->fatal("EPState::getTensorVar Out of Tensortial Var's range %d!", getNTensorVar() );
	 exit(1);
      }


}

//================================================================================
// Return Scalar pointer 
//================================================================================
double * EPState::getScalarVar() {
      
      return ScalarVar; 
      
}

//================================================================================
double * EPState::getScalarVar_commit() {
      
      return ScalarVar_commit; 
      
}

//================================================================================
double * EPState::getScalarVar_init() {
      
      return ScalarVar_init; 
      
}

//================================================================================
// Return Tensor pointer 
//================================================================================
stresstensor * EPState::getTensorVar() { 
    
      return TensorVar; 

}

//================================================================================
stresstensor * EPState::getTensorVar_commit() { 
    
      return TensorVar_commit; 

}

//================================================================================
stresstensor * EPState::getTensorVar_init() { 
    
      return TensorVar_init; 

}

//================================================================================
// set nth Scalar Var.... Starting from 1!!
//================================================================================
void EPState::setScalarVar(int WhichOne, double rval) { 

      if ( WhichOne <= getNScalarVar() )  
         ScalarVar[ WhichOne - 1 ] = rval; 
      else 
      {
         g3ErrorHandler->fatal("EPState::setScalarVar Out of ScalarVar's range %d!", getNScalarVar() );
         //cout << " Out of ScalarVar's range!";	  
	 exit(1);
      }
      
      
}

//================================================================================
// set nth Tensor Var.... Starting from 1!!
//================================================================================
void EPState::setTensorVar(int WhichOne, const stresstensor &rval) { 

      if ( WhichOne <= getNTensorVar() )  
         TensorVar[ WhichOne - 1 ] = rval;
      else 
      {
         g3ErrorHandler->fatal("EPState::setTensorVar Out of Tensor Var's range %d!", getNTensorVar() );	  
	 exit(1);
      }

}
//================================================================================
// set all Scalar Vars ..No boundary checking!
//================================================================================
void EPState::setScalarVar(double *rval) {

      if ( !rval ) {
         g3ErrorHandler->fatal("EPState::setScalarVar No scalar vars supplied");
         ::exit(1);  
      }
      for (int i = 0; i < getNScalarVar(); i++) {
         //cout << Scalarp[i] << endln; 
	 ScalarVar[i] = rval[i];
      }

}

//================================================================================
// set all Scalar Vars
//================================================================================
void EPState::setTensorVar(const stresstensor *rval) {

      if ( !rval ) {
         g3ErrorHandler->fatal("EPState::setTensorVar No tensorial vars supplied");
         ::exit(1);  
      }
      for (int i = 0; i < getNTensorVar(); i++) {
	 TensorVar[i] = rval[i];
 	 TensorVar[i].null_indices();
      }

}

//================================================================================
void EPState::print() { 
      cout << *this;

}


//================================================================================
// Set all state variables to initials

void EPState::setInit() { 
      
      stresstensor Stress0;
      straintensor Strain0;

      CurrentStress   = Stress_init;
      CurrentStrain   = Strain_init;
      ElasticStrain   = Strain0;
      PlasticStrain   = Strain0;
      dElasticStrain  = Strain0;
      dPlasticStrain  = Strain0;
      Eep = Eep_init;

      Stress_commit   = Stress_init;
      Strain_commit   = Strain_init;
      Eep_commit = Eep_init;

      for (int i = 0; i < NScalarVar; i++) {
          ScalarVar[i] = ScalarVar_init[i];
          ScalarVar_commit[i] = ScalarVar_init[i];
      }

      for (int i = 0; i < NTensorVar; i++) {
      	 TensorVar[i] = TensorVar_init[i];
      	 TensorVar_commit[i] = TensorVar_init[i];
      }

      Converged = false;	      	      			       
}

//================================================================================
int EPState::commitState () {

      int err = 0;
      // commit the variables state
      CurrentStress   = Stress_init;
      CurrentStrain   = Strain_init;
      Eep = Eep_init;

      Stress_commit   = CurrentStress;
      Strain_commit   = CurrentStrain;
      Eep_commit = Eep;

      for (int i = 0; i < NScalarVar; i++) {
          //ScalarVar[i] = ScalarVar_init[i];
          ScalarVar_commit[i] = ScalarVar[i];
      }

      for (int i = 0; i < NTensorVar; i++) {
      	 //TensorVar[i] = TensorVar_init[i];
      	 TensorVar_commit[i] = TensorVar[i];
      }

      return err;

}

//================================================================================
int EPState::revertToLastCommit () {
      int err = 0;
      // revert the variables state  to last commited
      CurrentStress   = Stress_commit;
      CurrentStrain   = Strain_commit;
      Eep = Eep_commit;
	     
      for (int i = 0; i < NScalarVar; i++) {
          //ScalarVar[i] = ScalarVar_init[i];
          ScalarVar[i] = ScalarVar_commit[i];
      }

      for (int i = 0; i < NTensorVar; i++) {
      	 //TensorVar[i] = TensorVar_init[i];
      	 TensorVar[i] = TensorVar_commit[i];
      }

      return err;

}

//================================================================================
int EPState::revertToStart () {

      int err = 0;
      // revert the variables state to the initials
      CurrentStress   = Stress_init;
      CurrentStrain   = Strain_init;
      Eep = Eep_init;

      Stress_commit   = Stress_init;
      Strain_commit   = Strain_init;
      Eep_commit = Eep_init;

      for (int i = 0; i < NScalarVar; i++) {
          ScalarVar[i] = ScalarVar_init[i];
          ScalarVar_commit[i] = ScalarVar_init[i];
      }

      for (int i = 0; i < NTensorVar; i++) {
      	 TensorVar[i] = TensorVar_init[i];
      	 TensorVar_commit[i] = TensorVar_init[i];
      }

      return err;
}





//================================================================================
// Overloaded Insertion Operator
// prints an EPState's contents 
//================================================================================
ostream & operator<< (ostream& os, const EPState & EPS)
    {
        os.setf( ios::showpos | ios::scientific);
        os.precision(4);
        os.width(10);       
        os << endln << "Elastic plastic state parameters: "  << endln;

        //os.width(10);       
        os << "\tEo = " << EPS.getEo() << ";";
        os << " E_Young = " << EPS.getE() << ";";
        //os.width(10);       
	os << " nu_Poisson = " << EPS.getnu() << ";";
	//os.width(10);       
	os << " rho = " << EPS.getrho() << endln;

        //os.width(10);       
        os << endln << "\tCurrent Stress:" << EPS.getStress() << endln;
        os << "\tIterati Stress:" << EPS.getIterativeStress() << endln;

	os << "\tCurrent Strain:" << EPS.getStrain() << endln;
	os << "\tElasticStrain :" << EPS.getElasticStrain() << endln;
	os << "\tPlasticStrain :" << EPS.getPlasticStrain() << endln;
	os << "\tdElasticStrain:" << EPS.getdElasticStrain() << endln;
	os << "\tdPlasticStrain:" << EPS.getdPlasticStrain() << endln;
	os << "\tEep.rank():" << EPS.getEep().rank() << endln;

        os.unsetf( ios::showpos );
	int NS = EPS.getNScalarVar();
	int NT = EPS.getNTensorVar();
	
	os << endln << "\tNScalarVar = " << NS << endln; 
    
        for (int i = 0; i < NS; i++) {
            os << "\tNo." << i+1 << " " << EPS.ScalarVar[i] << "; ";
	}
        os << endln << endln;
    
        os << "\tNTensorVar = " << NT;
        for (int i = 0; i < NT; i++) {
           os.unsetf( ios::showpos);
           os << endln << "\tNo." << i+1 << " tensorial var:";
           os.setf( ios::showpos);
           cout << EPS.TensorVar[i];
        }

        os << endln;           
        return os;
    }  

#endif

