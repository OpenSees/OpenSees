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
//# DATE:              19AUg2003
//# UPDATE HISTORY:    Sept 2003
//#		       May28, 2004
//#		       
//===============================================================================

#ifndef FiniteDeformationElastic3D_h
#define FiniteDeformationElastic3D_h

#include <math.h>

#include <ID.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <ConsoleErrorHandler.h>

#include <Matrix.h>
#include <Vector.h>
#include <Tensor.h>
#include <stresst.h>
#include <straint.h>

#include <NDMaterial.h>

class FiniteDeformationElastic3D : public NDMaterial
{
  public:
    
    FiniteDeformationElastic3D(int tag, int classTag, double );
    FiniteDeformationElastic3D();      
    
    virtual ~FiniteDeformationElastic3D();

    virtual double getRho(void);

    virtual int setTrialF(const straintensor &f);
    virtual int setTrialFIncr(const straintensor &df);
    virtual int setTrialC(const straintensor &c);
    virtual int setTrialCIncr(const straintensor &dc);

    virtual const Tensor& getTangentTensor(void) ;	  // Default Lagrangian Tangent Tensor
    virtual const Tensor& getInitialTangentTensor(void) ;

    virtual const  straintensor getStrainTensor(void) ;   // Default Green Lagrangian Strain
    virtual const  stresstensor getStressTensor(void) ;   // Default 2nd Piola Kirchhoff Stress
    virtual const  straintensor getF(void);
    virtual const  straintensor getC(void);

    virtual int commitState(void) ;
    virtual int revertToLastCommit(void) ;
    virtual int revertToStart(void) ;

    virtual NDMaterial *getCopy (void);
    virtual NDMaterial *getCopy (const char *type);

    virtual const char *getType (void) const;
    virtual int getOrder (void) const;

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    virtual void Print(OPS_Stream &s, int flag = 0);

    virtual int setParameter(char **argv, int argc, Information &info);
    virtual int updateParameter(int parameterID, Information &info);

    virtual const  stresstensor getPK1StressTensor(void) ;
    virtual const  stresstensor getCauchyStressTensor(void) ;

  protected:

    double rho;

};

#endif

