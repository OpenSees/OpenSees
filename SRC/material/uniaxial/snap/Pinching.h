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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:01:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/snap/Pinching.h,v $
//
//
// Pinching.h: implementation of the Pinching class from Fortran version.
// Originally from SNAP PROGRAM by Prof H.K. Krawinkler
//
// Written: A. Altoontash & Prof. G. Deierlein 12/01
// Revised: 03/02
//
// Purpose: This file contains the implementation for the Pinching class.
//
//////////////////////////////////////////////////////////////////////

// Pinching.h: interface for the Pinching class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PINCHING_H
#define PINCHING_H

#include <UniaxialMaterial.h>
#include <Vector.h>
#include <stdio.h>

class Pinching : public UniaxialMaterial  
{
 public:
  Pinching();
  Pinching(int tag, Vector inputParam );
  virtual ~Pinching();
  
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
  int recvSelf(int commitTag, Channel &theChannel,
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
 protected:
  void envelPosCap(double fy, double alfaPos, double alfaCap,
		   double cpDsp, double d, double *f, double *ek );
  
  void envelNegCap(double fy, double alfaNeg, double alfaCap,
		   double cpDsp, double d, double *f, double *ek);
  
  void recordInfo(int cond =0);
  
  
 private:
  
  // Input parameters
  double elstk,fyieldPos,fyieldNeg,alpha,Resfac; // Properties
  double capSlope,capDispPos,capDispNeg;	 // Cap
  double ecaps,ecapk,ecapa,ecapd,cs,ck,ca,cd;	 // Degradation parameters
  double fpPos,fpNeg,a_pinch;                    // Pinching
  
  // Parameters calculated from input data
  double ekhard,dyieldPos,dyieldNeg,Enrgts,Enrgta,Enrgtk,Enrgtd;
  double fPeakPos,fPeakNeg;
  double alfaNeg,alfaPos,fpdegPos,fpdegNeg,cpPos,cpNeg,fCapRefPos,fCapRefNeg;
  int flagDir;
  
  double hstv[18], hsLastCommit[18];
  
  FILE *OutputFile;		// For debugging
};

#endif
