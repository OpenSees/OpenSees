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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasHyster2Material.cpp,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasHyster2Material. FedeasHyster2Material wraps the FEDEAS
// 1d material subroutine Hyster_2.

#include <stdlib.h>
#include <FedeasHyster2Material.h>

FedeasHyster2Material::FedeasHyster2Material(int tag,
		double mom1p, double rot1p, double mom2p, double rot2p,
		double mom3p, double rot3p, double mom1n, double rot1n,
		double mom2n, double rot2n, double mom3n, double rot3n,
		double pinchX, double pinchY, double damfc1, double damfc2):
// 6 history variables and 16 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasHysteretic2, 6, 16)
{
	data[0]  = mom1p;
	data[1]  = rot1p;
	data[2]  = mom2p;
	data[3]  = rot2p;
	data[4]  = mom3p;
	data[5]  = rot3p;

	data[6]  = mom1n;
	data[7]  = rot1n;
	data[8]  = mom2n;
	data[9]  = rot2n;
	data[10] = mom3n;
	data[11] = rot3n;

	data[12] = pinchX;
	data[13] = pinchY;
	data[14] = damfc1;
	data[15] = damfc2;
}

FedeasHyster2Material::FedeasHyster2Material(int tag,
	double mom1p, double rot1p, double mom2p, double rot2p,
	double mom1n, double rot1n, double mom2n, double rot2n,
	double pinchX, double pinchY, double damfc1, double damfc2):
// 6 history variables and 16 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasHysteretic2, 6, 16)
{
	data[0]  = mom1p;
	data[1]  = rot1p;
	data[2]  = 0.5*(mom1p+mom2p);
	data[3]  = 0.5*(rot1p+rot2p);
	data[4]  = mom2p;
	data[5]  = rot2p;

	data[6]  = mom1n;
	data[7]  = rot1n;
	data[8]  = 0.5*(mom1n+mom2n);
	data[9]  = 0.5*(rot1n+rot2n);
	data[10] = mom2n;
	data[11] = rot2n;

	data[12] = pinchX;
	data[13] = pinchY;
	data[14] = damfc1;
	data[15] = damfc2;
}

FedeasHyster2Material::FedeasHyster2Material(int tag, const Vector &d):
// 6 history variables and 16 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasHysteretic2, 6, 16)
{
  if (d.Size() != numData) {
    opserr << "FedeasHyster2Material::FedeasHyster2Material -- not enough input arguments\n";
    exit(-1);
  }

  for (int i = 0; i < numData; i++)
    data[i] = d(i);
}

FedeasHyster2Material::FedeasHyster2Material(void):
FedeasMaterial(0, MAT_TAG_FedeasHysteretic2, 6, 16)
{
	// Does nothing
}

FedeasHyster2Material::~FedeasHyster2Material(void)
{
	// Does nothing
}

UniaxialMaterial*
FedeasHyster2Material::getCopy(void)
{
  Vector d(data, numData);

  FedeasHyster2Material *theCopy = new FedeasHyster2Material(this->getTag(), d);
  
  // Copy history variables
  for (int i = 0; i < 2*numHstv; i++)
    theCopy->hstv[i] = hstv[i];
  
  theCopy->epsilonP = epsilonP;
  theCopy->sigmaP   = sigmaP;
  theCopy->tangentP = tangentP;
  
  return theCopy;
}

double
FedeasHyster2Material::getInitialTangent(void)
{
	//return mom1p/rot1p;
	return data[0]/data[1];
}
