/*
################################################################################
# COPYRIGHT (C):     :-))                                                      #
# PROJECT:           Object Oriented Finite Element Program                    #
# PURPOSE:           General platform for elaso-plastic constitutive model     #
#                    implementation                                            #
# CLASS:             Template3Dep (the base class for all material point)      #
#                                                                              #
# VERSION:                                                                     #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
# TARGET OS:         DOS || UNIX || . . .                                      #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                               #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                               #
#                                                                              #
#                                                                              #
# DATE:              08-03-2000                                                #
# UPDATE HISTORY:    09-12-2000                                                #
#       May 2004, Zhao Cheng spliting the elastic part        #
#                                                                              #
#                                                                              #
# SHORT EXPLANATION: The Template3Dep class is used to hold specific           #
#                    yield surface, potential surface, Evolution law(s)        #
#                    and EPState of a 3D elasto-plastic material model for one #
#                    gauss point!! It is worthwhile noting that one model may  #
#                    have multiple evolution law.  Each evlotuion law is       #
#                    used to evolve one internal var.                          #
#                                                                              #
#                                                                              #
################################################################################
*/

#ifndef Template3Dep_H
#define Template3Dep_H

#include <NDMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <stresst.h>
#include <straint.h>
#include <BJtensor.h>

//#include <MD_YS.h>
//#include <MD_PS.h>
#include <YS.h>
#include <PS.h>
#include <EL_S.h>
#include <EL_T.h>
#include <EPState.h>

#include <Channel.h>
#include <G3Globals.h>

//#include <DP_YS.h>
//#include <DP_PS.h>
//#include <EPState.h>

//** Include the Elastic Material Models here
#include <ElasticIsotropic3D.h>
#include <ElasticCrossAnisotropic.h>
#include <PressureDependentElastic3D.h>

//#include <CDriver.h>

//class YieldSurface;
#define ITMAX 30
#define MAX_STEP_COUNT 30
#define NUM_OF_SUB_INCR 30
#define KK 1000.0  

class Template3Dep : public NDMaterial
{
 public:

  public:
    // constructor
    Template3Dep( int tag                ,
                  NDMaterial    &theElMat,
                  YieldSurface     *YS_ ,        
                  PotentialSurface *PS_ ,
                 EPState          *EPS_,
           EvolutionLaw_S   *ELS1_ , 
           EvolutionLaw_S   *ELS2_ , 
           EvolutionLaw_S   *ELS3_ , 
           EvolutionLaw_S   *ELS4_ , 
           EvolutionLaw_T   *ELT1_ ,
           EvolutionLaw_T   *ELT2_ ,
           EvolutionLaw_T   *ELT3_ ,
           EvolutionLaw_T   *ELT4_  );
    
    // Constructor0
    // If no evolution law is provided, then there will be no hardening or softening!
    Template3Dep(  int tag               ,
                   NDMaterial     &theElMat,
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
                  EPState          *EPS_);

    // Constructor1
    // Only one scalar evolution law is provided!
    Template3Dep(  int tag               ,
                   NDMaterial     &theElMat,
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
                  EPState          *EPS_,
            EvolutionLaw_S   *ELS1_ );

    // Constructor2
    // Only one tensorial evolution law is provided!
    Template3Dep(  int tag               ,
                   NDMaterial     &theElMat,
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
                  EPState          *EPS_,
            EvolutionLaw_T   *ELT1_ );

    // Constructor 3
    // One scalar evolution law and one tensorial evolution law are provided!
    Template3Dep(  int tag               , 
                   NDMaterial     &theElMat,
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
                  EPState          *EPS_,
            EvolutionLaw_S   *ELS1_, 
            EvolutionLaw_T   *ELT1_ );
    
    // Constructor 4
    // Two scalar evolution laws and one tensorial evolution law are provided!
    Template3Dep(  int tag               ,
                   NDMaterial     &theElMat,
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
                  EPState          *EPS_,
            EvolutionLaw_S   *ELS1_, 
            EvolutionLaw_S   *ELS2_, 
            EvolutionLaw_T   *ELT1_ );

    // Constructor 5
    // Two scalar evolution laws and two tensorial evolution laws are provided!
    Template3Dep(  int tag               ,
                   NDMaterial     &theElMat,
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
                  EPState          *EPS_,
            EvolutionLaw_S   *ELS1_, 
            EvolutionLaw_S   *ELS2_, 
            EvolutionLaw_T   *ELT1_,
            EvolutionLaw_T   *ELT2_ );

    // For parallel processing
    Template3Dep(void);
    virtual ~Template3Dep(void);

    // methods to set state and retrieve state using Matrix and Vector classes
    int setTrialStrain(const Vector &v);
    int setTrialStrain(const Vector &v, const Vector &r);
    int setTrialStrainIncr(const Vector &v) ;
    int setTrialStrainIncr(const Vector &v, const Vector &r) ;
    const Matrix &getTangent(void) ;
    const Matrix &getInitialTangent(void) ;

    const Vector &getStress(void) ;
    const Vector &getStrain(void) ;

    // methods to set and retrieve state using the Tensor class    
    int setTrialStrain(const Tensor &v) ;
    int setTrialStrain(const Tensor &v, const Tensor &r) ;    
    int setTrialStrainIncr(const Tensor &v) ;
    int setTrialStrainIncr(const Tensor &v, const Tensor &r) ;
    const Tensor &getTangentTensor(void) ;
    const stresstensor getStressTensor(void) ;
    const straintensor getStrainTensor(void) ;
    const straintensor getPlasticStrainTensor(void); //Added Joey Aug. 13, 2001
    double getpsi(void); //Added Joey 02-18-03

    EPState * getEPS() const;
    void setEPS( EPState &eps);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    
    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *code) ;

    //Template3Dep getCopy(void);  //????/
    //Template3Dep getCopy(const char *code) ;///???/

    const char *getType(void) const ;
    int getOrder(void) const ;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

    void Print(OPS_Stream &s, int flag =0);

    //Private Utility method
  //private:
    
     //These are from formerly CDriver
     EPState ForwardEulerEPState( const straintensor &strain_increment);
     
     EPState SemiBackwardEulerEPState( const straintensor &strain_increment);
     
     EPState FESubIncrementation( const straintensor &strain_increment,
                                  int number_of_subincrements);

     EPState BackwardEulerEPState( const straintensor &strain_increment);

     EPState BESubIncrementation( const straintensor & strain_increment,
                                  int number_of_subincrements);                                                 
  private:
                                                        
    //================================================================================
    // this one is intended to shell the previous three and to decide 
    // ( according to the data stored in Material_Model object ) 
    // which constitutive tensor to return ( forward ( non-constistent
    // or backward ( consistent ) or . . . 
    
    //virtual tensor ConstitutiveTensor(stresstensor   & final_stress, 
    //                                 stresstensor   & start_stress,
    //                                 straintensor   & strain_increment,
    //                                 Material_Model & Criterion,
    //                                 double           just_this_PP );
    

    //================================================================================
    // trying to find intersection point
    // according to M. Crisfield's book
    // "Non-linear Finite Element Analysis of Solids and Structures "
    // chapter 6.6.1 page 168.
    //================================================================================
    EPState PredictorEPState(straintensor & strain_increment);

    stresstensor yield_surface_cross(const stresstensor & start_stress,
                                     const stresstensor & end_stress);

    double zbrentstress(const stresstensor & start_stress,
                        const stresstensor & end_stress,
                        double x1, double x2, double tol);

    double func( const stresstensor &start_stress,
                 const stresstensor &end_stress,
                 double alfa );
    
   public: 
    tensor ElasticComplianceTensor(void) const;
    tensor ElasticStiffnessTensor(void) const;

   private:
    NDMaterial * getElMat() const;
    YieldSurface * getYS() const;
    PotentialSurface * getPS() const;

    //EPState &getEPS(void);
    
    //get scalar evolution laws
    EvolutionLaw_S * getELS1() const;
    EvolutionLaw_S * getELS2() const;
    EvolutionLaw_S * getELS3() const;    
    EvolutionLaw_S * getELS4() const;

    //get tensorial evolution laws
    EvolutionLaw_T * getELT1() const;
    EvolutionLaw_T * getELT2() const;
    EvolutionLaw_T * getELT3() const;
    EvolutionLaw_T * getELT4() const;
    

    // Get n copies of the NDMaterial 
    //NDMaterial **getCopy(int n) ;

    //================================================================================
    // Overloaded Insertion Operator
    //================================================================================
    friend OPS_Stream& operator<< (OPS_Stream& os, const Template3Dep & MP);

   private:

    NDMaterial  *theElasticMat;

    YieldSurface *YS;

    PotentialSurface *PS;

    EPState *EPS;
                               
    //Scalar variable evolution laws (currently at most 4  allowed)
    EvolutionLaw_S *ELS1; 
    EvolutionLaw_S *ELS2; 
    EvolutionLaw_S *ELS3; 
    EvolutionLaw_S *ELS4; 
    
    //Tensorial variable evolution laws (currently at most 4  allowed)
    EvolutionLaw_T *ELT1; 
    EvolutionLaw_T *ELT2; 
    EvolutionLaw_T *ELT3; 
    EvolutionLaw_T *ELT4; 

};


#endif

