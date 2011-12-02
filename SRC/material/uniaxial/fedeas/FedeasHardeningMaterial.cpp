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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-04-02 22:02:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasHardeningMaterial.cpp,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasHardeningMaterial. FedeasHardeningMaterial wraps the FEDEAS
// 1d material subroutine Hard_1.

#include <stdlib.h>
#include <FedeasHardeningMaterial.h>

FedeasHardeningMaterial::FedeasHardeningMaterial(int tag,
					 double E, double sigmaY, double Hiso, double Hkin):
// 3 history variables and 4 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasHardening, 3, 4)
{
	// Fill in material parameters
	data[0] = E;
	data[1] = sigmaY;
	data[2] = Hiso;
	data[3] = Hkin;
}

FedeasHardeningMaterial::FedeasHardeningMaterial(int tag, const Vector &d):
// 3 history variables and 4 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasHardening, 3, 4)
{
  if (d.Size() != numData) {
    opserr << "FedeasHardeningMaterial::FedeasHardeningMaterial -- not enough input arguments\n";
    exit(-1);	
  }

  for (int i = 0; i < numData; i++)
    data[i] = d(i);
}

FedeasHardeningMaterial::FedeasHardeningMaterial(void):
FedeasMaterial(0, MAT_TAG_FedeasHardening, 3, 4)
{
	// Does nothing
}

FedeasHardeningMaterial::~FedeasHardeningMaterial(void)
{
	// Does nothing
}

UniaxialMaterial*
FedeasHardeningMaterial::getCopy(void)
{
  Vector d(data, numData);

  FedeasHardeningMaterial *theCopy = new FedeasHardeningMaterial(this->getTag(), d);
  
  // Copy history variables
  for (int i = 0; i < 2*numHstv; i++)
    theCopy->hstv[i] = hstv[i];
  
  theCopy->epsilonP = epsilonP;
  theCopy->sigmaP   = sigmaP;
  theCopy->tangentP = tangentP;
  
  return theCopy;
}

double
FedeasHardeningMaterial::getInitialTangent(void)
{
	//return E;
	return data[0];
}
