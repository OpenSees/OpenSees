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
// $Date: 2007-10-26 04:32:06 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PrestressedSteelMaterial.h,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// PrestressedSteelMaterial. 

#ifndef PrestressedSteelMaterial_h
#define PrestressedSteelMaterial_h

#include <UniaxialMaterial.h>

#define MAT_TAG_PrestressedSteelMaterial 1976

class PrestressedSteelMaterial : public UniaxialMaterial
{
 public:
  PrestressedSteelMaterial(int tag, double E, double fy, double epslp, double epsd);
  PrestressedSteelMaterial();    
  ~PrestressedSteelMaterial();

  const char *getClassType(void) const {return "PrestressedSteelMaterial";};

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

  double E;
  double fy;
  double elp;
  double ed;
};

#endif

