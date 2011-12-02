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

#ifndef NeoHookeanCompressible3D_h
#define NeoHookeanCompressible3D_h

#include <FiniteDeformationElastic3D.h>

class NeoHookeanCompressible3D : public FiniteDeformationElastic3D
{
  public:
    NeoHookeanCompressible3D(int tag, int classTag, double, double, double );
    NeoHookeanCompressible3D(int tag, double, double, double );
    NeoHookeanCompressible3D();    
    
    virtual ~NeoHookeanCompressible3D();

    double getRho(void);

    int setTrialF(const straintensor &f);
    int setTrialFIncr(const straintensor &df);
    int setTrialC(const straintensor &c);
    int setTrialCIncr(const straintensor &dc);

    const Tensor& getTangentTensor(void) ;	  // Default Lagrangian Tangent Tensor
    const Tensor& getInitialTangentTensor(void) ;

    const  straintensor getStrainTensor(void) ;   // Default Green Lagrangian Strain
    const  stresstensor getStressTensor(void) ;   // Default 2nd Piola Kirchhoff Stress
    const  straintensor getF(void);
    const  straintensor getC(void);

    int commitState(void) ;
    int revertToLastCommit(void) ;
    int revertToStart(void) ;

    NDMaterial *getCopy (void);
    NDMaterial *getCopy (const char *type);

    const char *getType (void) const;
    int getOrder (void) const;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

//    int setParameter(char **argv, int argc, Information &info);
//    int updateParameter(int parameterID, Information &info);

    const  stresstensor getPK1StressTensor(void) ;
    const  stresstensor getCauchyStressTensor(void) ;

   
  private:

    int ComputeTrials(void);
     
  private:   
     
    double rho;
    double K;
    double G;

    straintensor F;
    straintensor C;
    double J;
    straintensor Cinv;

    int FromForC;

    Tensor Stiffness;
    straintensor thisGreenStrain;
    stresstensor thisPK2Stress;

};

#endif

