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
                                                                        
// $Revision: 1.5 $
// $Date: 2004-07-15 21:36:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasConcr1Material.cpp,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasConcr1Material. FedeasConcr1Material wraps the FEDEAS
// 1d material subroutine Concr_1.

#include <stdlib.h>
#include <FedeasConcr1Material.h>

FedeasConcr1Material::FedeasConcr1Material(int tag,
					 double fc, double ec, double fu, double eu):
// 3 history variables and 4 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasConcrete1, 3, 4)
{
	data[0]  = fc;
	data[1]  = ec;
	data[2]  = fu;
	data[3]  = eu;

	tangent = 2.0*data[0]/data[1];
	tangentP = tangent;
}

FedeasConcr1Material::FedeasConcr1Material(int tag, const Vector &d):
// 3 history variables and 4 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasConcrete1, 3, 4)
{
  if (d.Size() != numData) {
    opserr << "FedeasConcr1Material::FedeasConcr1Material -- not enough input arguments\n";
    exit(-1);
  }

  for (int i = 0; i < numData; i++)
    data[i] = d(i);

  tangent = 2.0*data[0]/data[1];
  tangentP = tangent;
}

FedeasConcr1Material::FedeasConcr1Material(void):
FedeasMaterial(0, MAT_TAG_FedeasConcrete1, 3, 4)
{
	// Does nothing
}

FedeasConcr1Material::~FedeasConcr1Material(void)
{
	// Does nothing
}

UniaxialMaterial*
FedeasConcr1Material::getCopy(void)
{
  Vector d(data, numData);

  FedeasConcr1Material *theCopy = new FedeasConcr1Material(this->getTag(), d);
  
  // Copy history variables
  for (int i = 0; i < 2*numHstv; i++)
    theCopy->hstv[i] = hstv[i];
  
  theCopy->epsilonP = epsilonP;
  theCopy->sigmaP   = sigmaP;
  theCopy->tangentP = tangentP;

  theCopy->epsilon = epsilonP;
  theCopy->sigma = sigmaP;
  theCopy->tangent = tangentP;  
  
  return theCopy;
}

double
FedeasConcr1Material::getInitialTangent(void)
{
	//return 2.0*fc/ec;
	return 2.0*data[0]/data[1];
}
