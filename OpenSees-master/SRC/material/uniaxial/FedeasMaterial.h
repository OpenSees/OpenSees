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
                                                                        
// $Revision: 1.7 $
// $Date: 2004-07-15 21:34:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FedeasMaterial.h,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasMaterial. FedeasMaterial provides a FORTRAN interface
// for programming uniaxial material models, using the subroutine
// interface from the FEDEAS ML1D library, developed by F.C. Filippou.
//
// For more information visit the FEDEAS web page:
//    http://www.ce.berkeley.edu/~filippou/Research/fedeas.htm

#ifndef FedeasMaterial_h
#define FedeasMaterial_h

#include <UniaxialMaterial.h>

class FedeasMaterial : public UniaxialMaterial
{
 public:
  FedeasMaterial(int tag, int classTag, int numHV, int numData);
  virtual ~FedeasMaterial();
  
  virtual int setTrialStrain(double strain, double strainRate = 0.0);
  virtual int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0);
  virtual double getStrain(void);
  virtual double getStress(void);
  virtual double getTangent(void);
  virtual double getInitialTangent(void) = 0;
  
  virtual int commitState(void);
  virtual int revertToLastCommit(void);    
  virtual int revertToStart(void);        
  
  virtual UniaxialMaterial *getCopy(void) = 0;
  
  virtual int sendSelf(int commitTag, Channel &theChannel);  
  virtual int recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker);    
  
  virtual void Print(OPS_Stream &s, int flag = 0);
  
 protected:
  // Invokes the FORTRAN subroutine
  virtual int invokeSubroutine(int ist);
  
  double *data;		// Material parameters array
  double *hstv;		// History array: first half is committed, second half is trial
  
  int numData;		// Number of material parameters
  int numHstv;		// Number of history variables
  
  double epsilonP;	// Committed strain
  double sigmaP;	// Committed stress
  double tangentP;	// Committed tangent

  double epsilon;	// Trial strain
  double sigma;		// Trial stress
  double tangent;	// Trial tangent

 private:
};

#endif
