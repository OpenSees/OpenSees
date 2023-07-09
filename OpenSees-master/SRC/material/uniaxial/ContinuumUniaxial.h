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
                                                                        
// $Revision: 1.1 $
// $Date: 2007-10-26 04:29:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ContinuumUniaxial.h,v $

// Written: MHS
// Created: June 2002
//
// Description: This file contains the class definition of ContinuumUniaxial.
// The ContinuumUniaxial class is a wrapper class that performs static
// condensation on a three-dimensional material model to give the 11
// uniaxial stress component.

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h> 
#include <UniaxialMaterial.h>

class NDMaterial;

class ContinuumUniaxial: public UniaxialMaterial {
 public:
  ContinuumUniaxial(int tag, NDMaterial &theMat);
  ContinuumUniaxial(void);
  virtual ~ContinuumUniaxial(void);
  
  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void);
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void);
  
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  
  UniaxialMaterial *getCopy(void);
  
  void Print(OPS_Stream &s, int flag);
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel,
	       FEM_ObjectBroker &theBroker);
  
  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int setParameter (const char **argv, int argc, Parameter &param);
  double getStressSensitivity(int gradIndex, bool conditional);
  int commitSensitivity(double strainGradient, int gradIndex, int numGrads);
  // AddingSensitivity:END ///////////////////////////////////////////
  

 private:
  double strain11;
  
  // Trial strain values
  double Tstrain22;
  double Tstrain33;
  double Tgamma12;
  double Tgamma23;
  double Tgamma31;

  // Committed strain values
  double Cstrain22;
  double Cstrain33;
  double Cgamma12;
  double Cgamma23;
  double Cgamma31;
  
  double initialTangent;

  NDMaterial *theMaterial;
};
