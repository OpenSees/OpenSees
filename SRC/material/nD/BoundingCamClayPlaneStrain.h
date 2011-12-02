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
                                                                        
// Written: Chris McGann
//          February 2011

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

#include <BoundingCamClay.h>

class BoundingCamClayPlaneStrain : public BoundingCamClay {

  public : 

  //null constructor
  BoundingCamClayPlaneStrain();

  //full constructor
  BoundingCamClayPlaneStrain(int tag, double c, double bulk, double OCR, double mu_o,
							 double alpha, double lambda, double h, double m, double mDen);


  //destructor
  ~BoundingCamClayPlaneStrain();

  NDMaterial* getCopy();
  const char* getType() const;
  int getOrder() const;

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

  //this is mike's problem
  int setTrialStrain(const Tensor &v);
  int setTrialStrain(const Tensor &v, const Tensor &r);    
  int setTrialStrainIncr(const Tensor &v);
  int setTrialStrainIncr(const Tensor &v, const Tensor &r);

  private :

  // static vectors and matrices
  static Vector strain;
  static Vector stress;
  static Matrix tangent;


};


