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
//# PURPOSE:           Hyper-spherical Constraint 
//# CLASS:             HSConstraint
//#
//# VERSION:           0.61803398874989 (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Ritu Jain, Boris Jeremic
//# PROGRAMMER(S):     Ritu, Boris Jeremic
//#
//#
//# DATE:              14Mar2003
//# UPDATE HISTORY:   
//#
//#
//===============================================================================
#ifndef HSConstraint_h
#define HSConstraint_h

#include <StaticIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class Matrix;

class HSConstraint : public StaticIntegrator
{
  public:
    HSConstraint(double arcLength, double psi_u=1.0, double psi_f=1.0, double u_ref=1.0);

    ~HSConstraint();

    int newStep(void);    
    int update(const Vector &deltaU);
    int domainChanged(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
    
  private:
    double arcLength2;
    double psi_u2;
    double psi_f2;
    double u_ref2;
    Matrix *scalingMatrix;
    Vector *deltaUhat;
    Vector *deltaUbar;
    Vector *deltaU;
    Vector *deltaUstep;
    Vector *phat; // the reference load vector
    double deltaLambdaStep;
    double currentLambda;
    int signLastDeltaLambdaStep;
};

#endif

