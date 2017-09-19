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

// Written: MHS
// Created: Dec 2012
//
// Description: This file contains the class definition for 
// InitStressNDMaterial.  InitStressNDMaterial wraps a UniaxialMaterial
// and imposes an initial strain


#ifndef InitStressNDMaterial_h
#define InitStressNDMaterial_h

#include <NDMaterial.h>

class InitStressNDMaterial : public NDMaterial
{
  public:
    InitStressNDMaterial(int tag, NDMaterial &material, const Vector &sigInit, int ndim = 3); 
    InitStressNDMaterial();
    ~InitStressNDMaterial();
    
    const char *getClassType(void) const {return "InitStressNDMaterial";};

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

    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *type);
    const char *getType(void) const;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
    int setParameter(const char **argv, int argc, Parameter &param);

    // send mass density to element in dynamic analysis
	double getRho(void);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    const Vector &getStressSensitivity(int gradIndex, bool conditional);
    int commitSensitivity(const Vector &depsdh, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    NDMaterial *theMaterial;
    Vector epsInit;
    Vector sigInit;
    int numDim;
};


#endif

