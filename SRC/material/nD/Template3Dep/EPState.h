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
#		     May 2004, Zhao Cheng spliting the elastic part	                        #
#       MAR2005 Guanzhou adding support for  const straintensor & ElasticStrain_commitp
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

#include <OPS_Globals.h>

// Constants
#define MaxNScalarVar 4
#define MaxNTensorVar 4
#define FALSE 0
#define TRUE 1

class EPState
{
  public:
//ZC05/2004    // Elastic parameters
//ZC05/2004    double Eo;                    // Young's modulus when p = p_atmosphere  -- [in-plane] in terms of cross-anisotropic material
//ZC05/2004    double E_Young;            // Young's modulus -- [in-plane] in terms of cross-anisotropic material
//ZC05/2004    double nu_Poisson;            // Poisson's ratio -- [in-plane] in terms of cross-anisotropic materi          
//ZC05/2004    double rho_mass_density;      // Mass density               
            
    // Trial state
    stresstensor CurrentStress;   // Current trial stress  --total          
    straintensor CurrentStrain;    // Current trial strain  --total       
    stresstensor IterativeStress; // Iterative Stress           
    straintensor ElasticStrain;    // Current cumulative elastic strain      
    straintensor PlasticStrain;    // Current cumulative plastic strain      
    straintensor dElasticStrain;  // Current elastic strain increment       
    straintensor dPlasticStrain;  // Current plastic strain increment       
    tensor       Eep;             // Current Elastic plastic stifness tensor
    int NScalarVar;               //Actual Number of internal scalar vars 
    int NTensorVar;               //Actual Number of internal tensor vars 
    
    double       ScalarVar[ MaxNScalarVar ]; // scalar variable array for scalar hardening vars 
    //static stresstensor TensorVar[ MaxNTensorVar ]; // tensor variable array for tensor hardening vars 
    stresstensor TensorVar[ MaxNTensorVar ]; // tensor variable array for tensor hardening vars 
    //straintensor TensorVar[ MaxNTensorVar ]; // tensor variable array for tensor hardening vars 

    // Commited state
    stresstensor Stress_commit;   // Commited stress  --total             
    straintensor Strain_commit;   // Commited strain  --total             
    straintensor ElasticStrain_commit; //GZ Mar2005
    //straintensor PlasticStrain


    double       ScalarVar_commit[ MaxNScalarVar ]; // Commited scalar variable array for scalar hardening vars 
    //static stresstensor TensorVar_commit[ MaxNTensorVar ]; // Commited tensor variable array for tensor hardening vars     
    stresstensor TensorVar_commit[ MaxNTensorVar ]; // Commited tensor variable array for tensor hardening vars 
    tensor       Eep_commit;      // Current Elastic plastic stifness tensor

    //Initial state
    stresstensor Stress_init;     // Initial stress  --total
    straintensor Strain_init;     // Initial strain  --total
    
    double       ScalarVar_init[ MaxNScalarVar ]; // initial scalar variable array for scalar hardening vars 
    //static stresstensor TensorVar_init[ MaxNTensorVar ]; // initial tensor variable array for tensor hardening vars     
    stresstensor TensorVar_init[ MaxNTensorVar ]; // initial tensor variable array for tensor hardening vars 
    tensor       Eep_init;             // initial Elastic plastic stifness tensor

    bool Converged;      // Bool to indicate whether this is the converged EPState by current CDriver
    
    int integratorFlag; //Guanzhou Mar2005

//ZC05/2004    // Flag to indicate if elastic portion is pressure dependent isotropic, pressure independent isotropic, pressure 
//ZC05/2004    // independent cross-anisotropic or pressure dependentcross-anisotropic 
//ZC05/2004    // 1 == pressure dependent isotropic (default case, for soil)
//ZC05/2004    // 2 == pressure independent isotropic
//ZC05/2004    // 3 == pressure independent cross-anisotropic
//ZC05/2004    // 4 == pressure dependent cross-anisotropic
//ZC05/2004
//ZC05/2004    int Elasticflag;
//ZC05/2004
//ZC05/2004    //Parameters for elastic cross-anistropic material
//ZC05/2004    double Ev;            // Ev: Young's modulus in a vertical direction -- [out-of-plane]
//ZC05/2004    double nuhv;   // nuhv: Poisson's ratio for strain in the vertical direction due to a horizontal direct stress -- [out-of-plane]
//ZC05/2004    double Ghv;          // Ghv: Modulus of shear deformation in a vertical plane -- [out-of-plane]
//ZC05/2004
//ZC05/2004    double eo;           // Initial void ratio Joey 02-11-03
//ZC05/2004    double ec;           // Void ratio at critical state at po Joey 02-11-03
//ZC05/2004    double Lambda;       // Slope of critical state line  Joey 02-11-03
//ZC05/2004    double po;           // Reference pressure (100 kPa by default)  Joey 02-11-03
    double e;            // Void ratio  Joey 02-11-03
    double psi;          // State parameter  Joey 02-18-03
//ZC05/2004    double a;            // Exponent in E = Eo (p/p_atm)^a for nonlinear elastic model Joey 02-11-03
    
  public:
    //Normal Constructor--no parameters
    EPState();
    ~EPState();
    
//ZC05/2004     //Normal Constructor1
//ZC05/2004     EPState(double              Eod,
//ZC05/2004             double              Ed,
//ZC05/2004             double              nu,
//ZC05/2004              double              rho,
//ZC05/2004             const stresstensor &stressp,
//ZC05/2004             const straintensor &strainp,
//ZC05/2004             const straintensor &Estrainp,
//ZC05/2004             const straintensor &Pstrainp,
//ZC05/2004             const straintensor &dEstrainp,
//ZC05/2004             const straintensor &dPstrainp,
//ZC05/2004             int                 NScalarp,
//ZC05/2004              const double       *Scalarp,
//ZC05/2004             int                 NTensorp,
//ZC05/2004              const stresstensor *Tensorp,
//ZC05/2004              const tensor       &Eepp,
//ZC05/2004           const stresstensor &Stress_commitp,
//ZC05/2004           const straintensor &Strain_commitp,    
//ZC05/2004           const double       *Scalar_commitp,
//ZC05/2004           const stresstensor *Tensor_commitp, 
//ZC05/2004           const tensor       &Eep_commitp,
//ZC05/2004           const stresstensor &Stress_initp,    
//ZC05/2004           const straintensor &Strain_initp,   
//ZC05/2004           const double       *Scalar_initp,
//ZC05/2004           const stresstensor *Tensor_initp,
//ZC05/2004           const tensor       &Eep_initp, 
//ZC05/2004           bool                Convergedp,
//ZC05/2004       int                 Elasticflagp = 0,
//ZC05/2004       double              Evp   = 0.0,
//ZC05/2004       double              nuhvp = 0.0,
//ZC05/2004       double              Ghvp = 0.0,
//ZC05/2004       double              eop = 0.85, 
//ZC05/2004       double              ecp = 0.80, 
//ZC05/2004       double              Lam = 0.025, 
//ZC05/2004       double              pop = 100.0,
//ZC05/2004       double              ep  = 0.85, 
//ZC05/2004       double              psip = 0.05, 
//ZC05/2004       double              ap  = 0.5
//ZC05/2004       );
//ZC05/2004 
//ZC05/2004 
//ZC05/2004     //Normal Constructor11
//ZC05/2004     EPState(double              Eod,
//ZC05/2004             double              Ed,
//ZC05/2004             double              nu,
//ZC05/2004             double              rho,
//ZC05/2004             const stresstensor  stressp,       
//ZC05/2004             const straintensor  strainp, 
//ZC05/2004             const straintensor  Estrainp,
//ZC05/2004             const straintensor  Pstrainp,
//ZC05/2004       int                 NScalarp,
//ZC05/2004       const double       *Scalarp,
//ZC05/2004       int                 NTensorp,
//ZC05/2004       const stresstensor *Tensorp, 
//ZC05/2004       int                 Elasticflagp = 0,
//ZC05/2004       double     Evp   = 0.0,
//ZC05/2004       double              nuhvp = 0.0,
//ZC05/2004       double              Ghvp = 0.0,
//ZC05/2004       double              eop = 0.85,
//ZC05/2004       double              ecp = 0.80,
//ZC05/2004       double              Lam = 0.025, 
//ZC05/2004       double              pop = 100.0,
//ZC05/2004       double              ap = 0.5
//ZC05/2004       );
//ZC05/2004 
//ZC05/2004     //Normal Constructor2
//ZC05/2004     EPState(double              Eod,    
//ZC05/2004             double              Ed,
//ZC05/2004             double              nu,
//ZC05/2004       double              rho,
//ZC05/2004             int                 NScalarp,
//ZC05/2004       const double       *Scalarp,
//ZC05/2004       int                 NTensorp,
//ZC05/2004       const stresstensor *Tensorp,
//ZC05/2004       int                 Elasticflagp = 0,
//ZC05/2004       double     Evp   = 0.0,
//ZC05/2004       double              nuhvp = 0.0,
//ZC05/2004       double              Ghvp = 0.0,
//ZC05/2004       double              eop = 0.85,
//ZC05/2004       double              ecp = 0.80,
//ZC05/2004       double              Lam = 0.025, 
//ZC05/2004       double              pop = 100.0,
//ZC05/2004       double              ap = 0.5
//ZC05/2004       );

    //Normal Constructor1
    EPState(const stresstensor &stressp,
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
            const straintensor & ElasticStrain_commitp,
	    const double       *Scalar_commitp,
            const stresstensor *Tensor_commitp, 
            const tensor       &Eep_commitp,
            const stresstensor &Stress_initp,    
            const straintensor &Strain_initp,   
            const double       *Scalar_initp,
            const stresstensor *Tensor_initp,
            const tensor       &Eep_initp, 
            bool                Convergedp, 
	    double              ep  = 0.85,
	    double              psip = 0.05,
	    int flag = 0 );

    //Normal Constructor11
    EPState(const stresstensor  stressp,       
            const straintensor  strainp, 
            const straintensor  Estrainp,
            const straintensor  Pstrainp,
            int                 NScalarp,
            const double       *Scalarp,
            int                 NTensorp,
            const stresstensor *Tensorp, 
	    double              ep = 0.85,
	    double              psip = 0.05, 
	    int flag = 0); //Guanzhou

    //Normal Constructor2
    EPState(int                 NScalarp,
            const double       *Scalarp,
            int                 NTensorp,
            const stresstensor *Tensorp );

    EPState *newObj();                 //create a clone of itself
    EPState( const EPState &rhs );     // Copy constructor    
    const EPState &operator=(const EPState &rhs );  //Overloading assignment

//ZC05/2004    int getElasticflag() const;

//ZC05/2004    double getE() const;
//ZC05/2004    double getEo() const;
//ZC05/2004    double getEv() const;
//ZC05/2004    double getnu() const;
//ZC05/2004    double getnuhv() const;
//ZC05/2004    double getGhv() const;
//ZC05/2004    double getrho() const;

    int getNScalarVar() const;
    int getNTensorVar() const;
    bool getConverged() const;

//ZC05/2004    // Added Joey 02-12-03
//ZC05/2004    double geteo() const;
//ZC05/2004    double getec() const;
//ZC05/2004    double getLam() const;
    double gete() const;
    double getpsi() const; //state parameter
    int getIntegratorFlag() const; //Guanzhou Mar 2005
//ZC05/2004    double getpo() const;
//ZC05/2004    double geta() const;

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
    straintensor getElasticStrain_commit() const;
    double * getScalarVar_commit();
    double getScalarVar_commit(int i);
    stresstensor * getTensorVar_commit();
    stresstensor getTensorVar_commit(int i);
    tensor getEep_commit() const;

    //Get initial state vars
    stresstensor getStress_init() const;
    straintensor getStrain_init() const;
    double * getScalarVar_init();
    double getScalarVar_init(int i);
    stresstensor * getTensorVar_init();
    stresstensor getTensorVar_init(int i);
    tensor getEep_init() const;

//ZC05/2004    void setElasticflag( int efd );
//ZC05/2004    void setEo( double Eod );
//ZC05/2004    void setE( double Ey );
//ZC05/2004    void setnu( double nud );

//ZC05/2004    void setEv( double Evd );
//ZC05/2004    void setGhv( double Ghvd );
//ZC05/2004    void setnuhv( double nud );

    void setStress( const stresstensor &newstress );
    void setIterativeStress( const stresstensor &newstress );
    void setStrain( const straintensor &newstrain );

    void setStress_commit( const stresstensor &newstress );
    void setStrain_commit( const straintensor &newstrain );

    void setStress_init( const stresstensor &newstress );
    void setStrain_init( const straintensor &newstrain );

    void setElasticStrain( const straintensor &newstrain );
    void setPlasticStrain( const straintensor &newstrain );
    void setdElasticStrain( const straintensor &newstrain );
    void setdPlasticStrain( const straintensor &newstrain );
    void setEep( const tensor &);
    void setConverged( bool b);

//ZC05/2004    // Added Joey 02-12-03
//ZC05/2004    void seteo( double eod );
//ZC05/2004    void setec( double ecd );
//ZC05/2004    void setLam( double Lamd );
    void sete( double ed );
    void setpsi( double psid );
//ZC05/2004    void setpo( double pod );
//ZC05/2004    void seta( double ad );

    // WhichOne starts from 1 and ends up to  NScalarVar
    double getScalarVar( int WhichOne) const;
    stresstensor getTensorVar(int WhichOne) const;

    // return the pointers
    double *getScalarVar();
    stresstensor *getTensorVar();

    // WhichOne starts from 1 and ends up to NTensorVar
    void setNScalarVar(int rval);
    void setScalarVar(int WhichOne, double rval);
    void setScalarVar_commit(int WhichOne, double rval);
    void setScalarVar_init(int WhichOne, double rval);

    void setNTensorVar(int rval);
    void setTensorVar(int WhichOne, const stresstensor &rval);
    void setTensorVar_commit(int WhichOne, const stresstensor &rval);
    void setTensorVar_init(int WhichOne, const stresstensor &rval);

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
    friend OPS_Stream & operator<< (OPS_Stream& os, const EPState & EPS);
    //friend ostream & operator<< (ostream& os, const EPState & EPS);

};


#endif

