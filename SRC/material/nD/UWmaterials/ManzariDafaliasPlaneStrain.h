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
// Description: This file contains the class definition for ManzariDafaliasPlaneStrain.

#ifndef ManzariDafaliasPS_h
#define ManzariDafaliasPS_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

#include <ManzariDafalias.h>

class ManzariDafaliasPlaneStrain : public ManzariDafalias
{
  public : 

  //full constructor
  ManzariDafaliasPlaneStrain(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c,
	double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd,double z_max, 
	double cz, double massDen, double TolF = 1e-7, double TolR = 1e-7, int jacoType = 1, int integrationScheme = 1);
  //null constructor
  ManzariDafaliasPlaneStrain();

  //destructor
  ~ManzariDafaliasPlaneStrain();

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

  //send back the tangent 
  const Matrix& getTangent();
  const Matrix& getInitialTangent();
  
  //send back the state parameters
  const Vector  getState();

  private :

  // static vectors and matrices
  static Vector strain;
  static Vector stress;
  static Matrix tangent;

};

#endif