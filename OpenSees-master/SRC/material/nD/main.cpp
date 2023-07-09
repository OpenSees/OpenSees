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
// $Date: 2005-12-14 23:49:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/main.cpp,v $
                                                                        
                                                                        
#include "Vector.h"
#include "ID.h"
#include "Matrix.h"

#include <ElasticIsotropicThreeDimensional.h>

#include <OPS_Globals.h>
#include <StandardStream.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

int main()
{
  //                                      tag     E      nu  rho
  ElasticIsotropicThreeDimensional material(1, 29000.0, 0.3, 0.0);

  Vector strain(6); // 11,22,33,12,23,31
  strain(0) = 0.002;

  material.setTrialStrain(strain);
  const Vector &stress = material.getStress();
  const Matrix &tangent = material.getTangent();
  opserr << stress << tangent << endln;
}



