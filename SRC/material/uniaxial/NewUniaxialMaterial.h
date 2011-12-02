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
                                                                        
// $Revision: 1.6 $
// $Date: 2006-08-15 00:41:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/NewUniaxialMaterial.h,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// NewUniaxialMaterial. 

#ifndef NewUniaxialMaterial_h
#define NewUniaxialMaterial_h

#include <UniaxialMaterial.h>

#define MAT_TAG_NewUniaxialMaterial 1976

class NewUniaxialMaterial : public UniaxialMaterial
{
 public:
  NewUniaxialMaterial(int tag);
  NewUniaxialMaterial();    
  ~NewUniaxialMaterial();

  const char *getClassType(void) const {return "NewUniaxialMaterial";};

  int setTrialStrain(double strain, double strainRate = 0.0); 
  double getStrain(void);
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void);
  
  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);        
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
 protected:
  
 private:
  double trialStrain;   // trial strain
  double trialStress;   // trial stress
  double trialTangent;  // trial tangent
  double commitStrain;   // commit strain
  double commitStress;   // commit stress
  double commitTangent;  // commit tangent
};

#endif

