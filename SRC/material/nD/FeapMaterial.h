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
// $Date: 2006-08-04 18:18:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/FeapMaterial.h,v $
                                                                      
// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// FeapMaterial. FeapMaterial wraps a Feap material subroutine.

#ifndef FeapMaterial_h
#define FeapMaterial_h

#include <NDMaterial.h>

class FeapMaterial : public NDMaterial
{
 public:
  FeapMaterial(int tag, int classTag, int numHV, int numData,
	       double rho = 0.0);
  FeapMaterial(int classTag);
  virtual ~FeapMaterial();

  virtual const char *getClassType(void) const {return "FeapMaterial";};
  
  virtual int setTrialStrain(const Vector &strain);
  virtual const Vector &getStrain(void);
  virtual const Vector &getStress(void);
  virtual const Matrix &getTangent(void);
  virtual double getRho(void);
  
  virtual int commitState(void);
  virtual int revertToLastCommit(void);    
  virtual int revertToStart(void);        
  
  virtual NDMaterial *getCopy(void);
  virtual NDMaterial *getCopy(const char *type);
  virtual const char *getType(void) const;
  virtual int getOrder(void) const;
  
  virtual int sendSelf(int commitTag, Channel &theChannel);  
  virtual int recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker);    
  
  virtual void Print(OPS_Stream &s, int flag = 0);
  
 protected:
  // Invokes the Feap subroutine
  virtual int invokeSubroutine(int isw);
  virtual int fillDArray(void);
  
  double *ud;	// Material parameters array
  double *hstv;	// History array: first half is committed, second half is trial
  
  static double d[];    // Feap material parameter array
  
  double rho;           // Material mass density
  
 private:
  int numHV;		// Number of history variables
  int numData;		// Number of material parameters
  
  double eps[6];        // Strain vector
  static double sig[6]; // Stress vector
  static double dd[36]; // Tangent matrix
  
  static Vector strain3;
  static Vector strain4;
  static Vector strain6;
  
  static Vector sigma3;
  static Vector sigma4;
  static Vector sigma6;
  
  static Matrix tangent3;
  static Matrix tangent4;
  static Matrix tangent6;
  
  enum Formulation{Unknown, ThreeDimensional, PlaneStrain, AxiSymmetric};
  int myFormulation;
  void setType(Formulation form);
};

#endif
