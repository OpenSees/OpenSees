//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Finite Deformation Hyper-Elastic classes
//# CLASS:
//#
//# VERSION:           0.6_(1803398874989) (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              July 2004
//# UPDATE HISTORY:
//#
//===============================================================================

#ifndef FDdecoupledElastic3D_H
#define FDdecoupledElastic3D_H

#include <FiniteDeformationElastic3D.h>

#include <W.h>  // for W Strain Energy Functions

class FDdecoupledElastic3D : public FiniteDeformationElastic3D
{
  public:
    FDdecoupledElastic3D(int tag, int classTag, WEnergy * , double );
    FDdecoupledElastic3D(int tag, WEnergy * , double );
    FDdecoupledElastic3D(int tag, WEnergy * );
    FDdecoupledElastic3D();       
    
    virtual ~FDdecoupledElastic3D();

    FDdecoupledElastic3D(FDdecoupledElastic3D &fde3d); 

    const char *getClassType(void) const {return "FDdecoupledElastic3D";};

    double getRho(void);

    int setTrialF(const straintensor &f);
    int setTrialFIncr(const straintensor &df);
    int setTrialC(const straintensor &c);
    int setTrialCIncr(const straintensor &dc);

    const Tensor& getTangentTensor(void) ;	  // Default Lagrangian Tangent Tensor
    const Tensor& getInitialTangentTensor(void) ;

    const  straintensor& getStrainTensor(void) ;   // Default Green Lagrangian Strain
    const  stresstensor& getStressTensor(void) ;   // Default 2nd Piola Kirchhoff Stress
    const  straintensor& getF(void);
    const  straintensor& getC(void);

//    virtual const Vector &getStress(void);
//    virtual const Vector &getStrain(void);

//    virtual const stresstensor getCommittedStress(void);
//    virtual const straintensor getCommittedStrain(void);

//    virtual const straintensor getPlasticStrainTensor(void);

    int commitState(void) ;
    int revertToLastCommit(void) ;
    int revertToStart(void) ;

    NDMaterial *getCopy (void);
    NDMaterial *getCopy (const char *type);

    const char *getType (void) const;
    //int getOrder (void) const;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

//    int setParameter(char **argv, int argc, Information &info);
//    int updateParameter(int parameterID, Information &info);

    const  stresstensor& getPK1StressTensor(void) ;
    const  stresstensor& getCauchyStressTensor(void) ;

  private:
  //public:
    WEnergy *getWEnergy(void);

    const double getJ(void) ;
    const Vector getlambda(void) ;
    const Vector getlambda_wave(void) ;

    const Vector wa(void) ;
    const Tensor Yab(void) ;
    const Tensor FDisoStiffness(void) ;
    const Tensor FDvolStiffness(void) ;


//    int setInitialTangentTensor(void);
    int ComputeTrials(void);
    int getCaseIndex(void);

  protected:

     WEnergy * W;

     double rho;

     straintensor F;
     straintensor C;
     double J;
     straintensor Cinv;
     double lambda1, lambda2, lambda3;
     double lambda_wave1, lambda_wave2, lambda_wave3;
     int caseIndex; 
     int FromForC;

     Tensor Stiffness;
     straintensor thisGreenStrain;
     stresstensor thisPK2Stress;

     static stresstensor static_FDE_stress;  // only for static reference return
};

#endif

