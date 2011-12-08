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
  ManzariDafaliasPlaneStrain(int tag, double Ko, double Go, double v, double b, double Patm,
								      double Ao, double ho, double Cm, double Me, double Mc,
									  double kBE, double kBC, double kDE, double kDC, double ecRef,
									  double lambda, double Pref, double m, double Fmax, double Cf,
									  double eo, double mDen);
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

  private :

  // static vectors and matrices
  static Vector strain;
  static Vector stress;
  static Matrix tangent;

};
