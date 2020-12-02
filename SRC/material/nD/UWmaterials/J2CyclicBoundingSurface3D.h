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
                                                                       
// Written: UW Computational Geomechanics Group
//          Pedro Arduino (*), ALborz Ghofrani(*)
//			(*)  University of Washington
//          March 2020
//
// Description: This file contains the class definition for J2CyclicBoundingSurface3D.
#ifndef J2CyclicBoundingSurface3D_h
#define J2CyclicBoundingSurface3D_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <J2CyclicBoundingSurface.h>

class J2CyclicBoundingSurface3D : public J2CyclicBoundingSurface
{
  public : 

  //full constructor
  J2CyclicBoundingSurface3D(int tag, double G, double K, double su, double rho, double h, double m, double h0, double chi, double beta);

  //null constructor
  J2CyclicBoundingSurface3D();

  //destructor
  ~J2CyclicBoundingSurface3D();

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

  private :

};


#endif
