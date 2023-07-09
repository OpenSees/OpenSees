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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasSteel2Material.cpp,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasSteel2Material. FedeasSteel2Material wraps the FEDEAS
// 1d material subroutine Steel_2.

#include <stdlib.h>
#include <FedeasSteel2Material.h>

FedeasSteel2Material::FedeasSteel2Material(int tag,
					 double fy, double E0, double b,
					 double R0, double cR1, double cR2,
					 double a1, double a2, double a3, double a4):
// 8 history variables and 10 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasSteel2, 8, 10)
{
	data[0]  = fy;
	data[1]  = E0;
	data[2]  = b;
	data[3]  = R0;
	data[4]  = cR1;
	data[5]  = cR2;
	data[6]  = a1;
	data[7]  = a2;
	data[8]  = a3;
	data[9]  = a4;

	tangent = E0;
	tangentP = E0;
}

FedeasSteel2Material::FedeasSteel2Material(int tag,
					 double fy, double E0, double b,
					 double R0, double cR1, double cR2):
// 8 history variables and 10 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasSteel2, 8, 10)
{
	data[0]  = fy;
	data[1]  = E0;
	data[2]  = b;
	data[3]  = R0;
	data[4]  = cR1;
	data[5]  = cR2;

	// Default values for no isotropic hardening
	data[6]  = 0.0;
	data[7]  = 1.0;
	data[8]  = 0.0;
	data[9]  = 1.0;

	tangent = E0;
	tangentP = E0;
}

FedeasSteel2Material::FedeasSteel2Material(int tag,
					 double fy, double E0, double b):
// 8 history variables and 10 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasSteel2, 8, 10)
{
	data[0]  = fy;
	data[1]  = E0;
	data[2]  = b;

	// Default values for elastic to hardening transitions
	data[3]  = 15.0;
	data[4]  = 0.925;
	data[5]  = 0.15;

	// Default values for no isotropic hardening
	data[6]  = 0.0;
	data[7]  = 1.0;
	data[8]  = 0.0;
	data[9]  = 1.0;

	tangent = E0;
	tangentP = E0;
}

FedeasSteel2Material::FedeasSteel2Material(int tag, const Vector &d):
// 8 history variables and 10 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasSteel2, 8, 10)
{
  if (d.Size() != numData) {
    opserr << "FedeasSteel2Material::FedeasSteel2Material -- not enough input arguments\n";
    exit(-1);
  }

  for (int i = 0; i < numData; i++)
    data[i] = d(i);
}

FedeasSteel2Material::FedeasSteel2Material(void):
FedeasMaterial(0, MAT_TAG_FedeasSteel2, 8, 10)
{
	// Does nothing
}

FedeasSteel2Material::~FedeasSteel2Material(void)
{
	// Does nothing
}

UniaxialMaterial*
FedeasSteel2Material::getCopy(void)
{
  Vector d(data, numData);

  FedeasSteel2Material *theCopy = new FedeasSteel2Material(this->getTag(), d);
  
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
FedeasSteel2Material::getInitialTangent(void)
{
	//return E;
	return data[1];
}
