/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision$
// $Date$
// $Source$

// Description: This file contains the class definition for 
// MinMaxNDMaterial.  MinMaxNDMaterial wraps an NDMaterial
// and imposes an initial strain
//
// MHS


#ifndef MinMaxNDMaterial_h
#define MinMaxNDMaterial_h

#include <NDMaterial.h>

class MinMaxNDMaterial : public NDMaterial
{
  public:
  MinMaxNDMaterial(int tag, NDMaterial &material, double epsmin, double epsmax);
    MinMaxNDMaterial();
    ~MinMaxNDMaterial();
    
    const char *getClassType(void) const { return "MinMaxNDMaterial"; }

    int setTrialStrain(const Vector &strain); 
    int setTrialStrain(const Vector &strain, const Vector &rate); 
    int setTrialStrainIncr(const Vector &strain); 
    int setTrialStrainIncr(const Vector &strain, const Vector &rate); 

    const Vector &getStress(void);
    const Vector &getStrain(void);          

    const Matrix &getTangent(void);
    const Matrix &getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    // send mass density to element in dynamic analysis
    double getRho(void);

    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *type);
    const char *getType(void) const;
    int getOrder(void) const;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information& info);

    Response* setResponse(const char** argv, int argc, OPS_Stream& output);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    const Vector &getStressSensitivity(int gradIndex, bool conditional);
    int commitSensitivity(const Vector &depsdh, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  private:
    NDMaterial *theMaterial;

  double minStrain;
  double maxStrain;

  bool Tfailed;
  bool Cfailed;

  enum type {PlaneStress2d = 1, BeamFiber = 2, BeamFiber2d = 3, BeamFiber2dPS = 4, PlateFiber = 5, ThreeDimensional = 6};
  type myType;
};


#endif

