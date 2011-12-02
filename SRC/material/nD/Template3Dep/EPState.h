/*
################################################################################
# COPYRIGHT (C):     :-))                                                      #
# PROJECT:           Object Oriented Finite Element Program                    #
# PURPOSE:           General platform for elaso-plastic constitutive model     #
#                    implementation                                            #
# CLASS:             DPEPState (the base class for all Elasto-plastic state)   #
#                                                                              #
# VERSION:                                                                     #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
# TARGET OS:         DOS || UNIX || . . .                                      #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                               #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                               #
#                                                                              #
#                                                                              #
# DATE:              08-03-2000                                                #
# UPDATE HISTORY:                                                              #
#                                                                              #
#                                                                              #
#                                                                              #
#                                                                              #
# SHORT EXPLANATION: This class is used to hold all state parameters and       #
#                    internal variables in an elasto-plastic constitutive      #
#                    model.                                                    #
################################################################################
*/

#ifndef EPState_H
#define EPState_H

#include <stresst.h>
#include <straint.h>
#include <BJtensor.h>

#include <iostream.h>
#include <iomanip.h>

#define endln "\n"

// Constants
#define MaxNScalarVar 4
#define MaxNTensorVar 4
#define FALSE 0
#define TRUE 1

class EPState
{
  // Elastic parameters
  public:
    double Eo;	                  // Young's modulus when p = p_atmosphere  
    double E_Young;	          // Young's modulus  			    
    double nu_Poisson;	          // Poisson's ratio  			    
    double rho_mass_density;      // Mass density     			    
    				
    // Trial state
    stresstensor CurrentStress;   // Current trial stress  --total          
    straintensor CurrentStrain;	  // Current trial strain  --total 	    
    stresstensor IterativeStress; // Iterative Stress 			    
    straintensor ElasticStrain;	  // Current cumulative elastic strain      
    straintensor PlasticStrain;	  // Current cumulative plastic strain      
    straintensor dElasticStrain;  // Current elastic strain increment       
    straintensor dPlasticStrain;  // Current plastic strain increment       
    tensor       Eep;             // Current Elastic plastic stifness tensor
    int NScalarVar;    	  	  //Actual Number of internal scalar vars 
    int NTensorVar;    	  	  //Actual Number of internal tensor vars 
    double       ScalarVar[ MaxNScalarVar ]; // scalar variable array for scalar hardening vars 
    stresstensor TensorVar[ MaxNTensorVar ]; // tensor variable array for tensor hardening vars 

    // Commited state
    stresstensor Stress_commit;   // Commited stress  --total               
    straintensor Strain_commit;	  // Commited strain  --total 	            
    double       ScalarVar_commit[ MaxNScalarVar ]; // Commited scalar variable array for scalar hardening vars 
    stresstensor TensorVar_commit[ MaxNTensorVar ]; // Commited tensor variable array for tensor hardening vars 
    tensor       Eep_commit;      // Current Elastic plastic stifness tensor

    //Initial state
    stresstensor Stress_init;     // Initial stress  --total               
    straintensor Strain_init;	  // Initial strain  --total 	            
    double       ScalarVar_init[ MaxNScalarVar ]; // initial scalar variable array for scalar hardening vars 
    stresstensor TensorVar_init[ MaxNTensorVar ]; // initial tensor variable array for tensor hardening vars 
    tensor       Eep_init;             // initial Elastic plastic stifness tensor

    bool Converged;      // Bool to indicate whether this is the converged EPState by current CDriver
    
  public:
    //Normal Constructor--no parameters
    EPState();                         
    
    //Normal Constructor1
    EPState(double              Eod,    
            double              Ed,
            double              nu,
	    double              rho,
            const stresstensor &stressp,
            const straintensor &strainp,
            const straintensor &Estrainp,
            const straintensor &Pstrainp,
            const straintensor &dEstrainp,
            const straintensor &dPstrainp,
            int                 NScalarp,
	    const double       *Scalarp,
            int                 NTensorp,
	    const stresstensor *Tensorp,
	    const tensor       &Eepp,
    	    const stresstensor &Stress_commitp,
    	    const straintensor &Strain_commitp,	  
    	    const double       *Scalar_commitp,
    	    const stresstensor *Tensor_commitp, 
    	    const tensor       &Eep_commitp,
    	    const stresstensor &Stress_initp,    
    	    const straintensor &Strain_initp,	 
    	    const double       *Scalar_initp,
    	    const stresstensor *Tensor_initp,
    	    const tensor       &Eep_initp, 
            bool                Convergedp      
	    );


    EPState(double             Eod,
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
	    const tensor     * Tensorp );

    //Normal Constructor2
    EPState(double              Eod,    
            double              Ed,
            double              nu,
	    double              rho,
            int                 NScalarp,
	    const double       *Scalarp,
	    int                 NTensorp,
	    const tensor       *Tensorp);

    EPState *newObj();                 //create a clone of itself
    EPState( const EPState &rhs );     // Copy constructor    
    const EPState & EPState::operator=(const EPState &rhs );  //Overloading assignment

    double getE() const;
    double getEo() const;
    double getnu() const;
    double getrho() const;
    int getNScalarVar() const;
    int getNTensorVar() const;
    bool getConverged() const;

    stresstensor getStress() const;
    stresstensor getIterativeStress() const;
    straintensor getStrain() const;
    straintensor getElasticStrain() const;
    straintensor getPlasticStrain() const;
    straintensor getdElasticStrain() const;
    straintensor getdPlasticStrain() const;
    tensor getEep() const;

    //Get commited state vars
    stresstensor getStress_commit() const;
    straintensor getStrain_commit() const;
    double *getScalarVar_commit();
    stresstensor *getTensorVar_commit();
    tensor getEep_commit() const;

    //Get initial state vars
    stresstensor getStress_init() const;
    straintensor getStrain_init() const;
    double *getScalarVar_init();
    stresstensor *getTensorVar_init();
    tensor getEep_init() const;

    void setE( double Ey );
    void setStress( const stresstensor &newstress );
    void setIterativeStress( const stresstensor &newstress );
    void setStrain( const straintensor &newstrain );
    void setElasticStrain( const straintensor &newstrain );
    void setPlasticStrain( const straintensor &newstrain );
    void setdElasticStrain( const straintensor &newstrain );
    void setdPlasticStrain( const straintensor &newstrain );
    void setEep( const tensor &);
    void setConverged( bool b);

    // WhichOne starts from 1 and ends up to  NScalarVar
    double getScalarVar( int WhichOne) const;
    stresstensor getTensorVar(int WhichOne) const;

    // return the pointers
    double *getScalarVar();
    stresstensor *getTensorVar();

    // WhichOne starts from 1 and ends up to NTensorVar
    void setScalarVar(int WhichOne, double rval);
    void setTensorVar(int WhichOne, const stresstensor &rval);
    void setScalarVar(double *rval);
    void setTensorVar(const stresstensor *rval);
    void setInit();

    //added for OpenSees _Zhaohui 02-10-2000
    int commitState ();
    int revertToLastCommit ();
    int revertToStart ();

    void print();

    //================================================================================
    // Overloaded Insertion Operator
    // prints an EPState's contents 
    //================================================================================
    friend ostream & operator<< (ostream& os, const EPState & EPS);

};


#endif

