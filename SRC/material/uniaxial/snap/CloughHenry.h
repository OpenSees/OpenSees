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
                                                                        
// CloughHenry.h: based on Clough.h/.cpp, slight modification
//
// Original Clough Written: Arash Altoontash,
// Modification: fmk
//
// Purpose: This file contains the implementation for the CloughHenry class.
//
//////////////////////////////////////////////////////////////////////

// CloughHenry.h: interface for the CloughHenry class.
//
//////////////////////////////////////////////////////////////////////

#ifndef CloughHenry_H
#define CloughHenry_H

#include <UniaxialMaterial.h>
#include <Vector.h>
#include <stdio.h>

class CloughHenry : public UniaxialMaterial
{
 public:
  CloughHenry();
  CloughHenry(int tag, Vector inputParam );
  virtual ~CloughHenry();
  const char *getClassType(void) const {return "CloughHenry";};	
  int setTrialStrain(double d, double strainRate = 0.0);
  double getStrain(void);
  
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void);
  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);  
	
  //virtual
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel,FEM_ObjectBroker &theBroker);
  
  void Print(OPS_Stream &s, int flag =0);
  
 protected:
  void envelPosCap(double fy, double alphaPos, double alphaCap,
		   double cpDsp, double d, double *f, double *ek );
  
  void envelNegCap(double fy, double alphaNeg, double alphaCap,
		   double cpDsp, double d, double *f, double *ek);
  
  void recordInfo(int cond =0);
  
  
 private:
  // Input parameters
  double elstk,fyieldPos,fyieldNeg,alpha,Resfac;		//	Properties
  double capSlope,capDispPos,capDispNeg;	 // Cap
  double ecaps,ecapk,ecapa,ecapd,cs,ck,ca,cd;	 // Degradation parameters
  
  // Parameters calculated from input data
  double dyieldPos,dyieldNeg;
  double Enrgts,Enrgtk,Enrgta,Enrgtd;
  
  double hsTrial[24], hsCommit[24], hsLastCommit[24];
  
};

#endif
