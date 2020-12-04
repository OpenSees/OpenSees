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
                                                                       
// Created: Pedro Arduino, UW, 11.2011
//
// Description: This file contains the class definition for ManzariDafaliasPlaneStrainRO.

#ifndef ManzariDafaliasPSRO_h
#define ManzariDafaliasPSRO_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

#include "ManzariDafaliasRO.h"

class ManzariDafaliasPlaneStrainRO : public ManzariDafaliasRO
{
  public : 

  //full constructor
  ManzariDafaliasPlaneStrainRO(int tag, double G0, double nu, double B, double a1, double gamma1, double e_init, double Mc, double c, 
					double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd, 
					double z_max, double cz, double mDen, double kappa = 2.0, int integrationScheme = 2, int tangentType = 2, int JacoType = 1, 
					double TolF = 1.0e-7, double TolR = 1.0e-7);
  //null constructor
  ManzariDafaliasPlaneStrainRO();

  //destructor
  ~ManzariDafaliasPlaneStrainRO();

  NDMaterial *getCopy(void);
  const char *getType(void) const;
  int getOrder(void) const;

  int setTrialStrain(const Vector &strain_from_element);

  // Unused trialStrain functions
  int setTrialStrain(const Vector &v, const Vector &r);
    
  //send back the strain
  const Vector& getStrain();

  //send back the stress 
  const Vector& getStress();
  const Vector& getStressToRecord();

  //send back the tangent 
  const Matrix& getTangent();
  const Matrix& getInitialTangent();
  
  //send back the state parameters
  const Vector  getState();

  private :

  // static vectors and matrices
  static Vector mSigma_M  ; // mSigma with continuum mechanic sign convention
  static Vector mEpsilon_M; // mEpsilon with continuum mechanic sign convention
  static Vector rSigma;     // Stress for the recorders
  static Matrix mTangent;
  static Matrix mTangent_init;

};

#endif
