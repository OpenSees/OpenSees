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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasConcr3Material.cpp,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasConcr3Material. FedeasConcr3Material wraps the FEDEAS
// 1d material subroutine Concr_3.

#include <stdlib.h>
#include <FedeasConcr3Material.h>

FedeasConcr3Material::FedeasConcr3Material(int tag,
					 double fc, double ec, double fu, double eu,
					 double rat, double ft, double epst0,
					 double ft0, double beta, double epstu):
// 2 history variables and 10 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasConcrete3, 2, 10)
{
	data[0]  = fc;
	data[1]  = ec;
	data[2]  = fu;
	data[3]  = eu;
	data[4]  = rat;
	data[5]  = ft;
	data[6]  = epst0;
	data[7]  = ft0;
	data[8]  = beta;
	data[9]  = epstu;
}

FedeasConcr3Material::FedeasConcr3Material(int tag, const Vector &d):
// 2 history variables and 10 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasConcrete3, 2, 10)
{
  if (d.Size() != numData) {
    opserr << "FedeasConcr3Material::FedeasConcr3Material -- not enough input arguments\n";
    exit(-1);
  }
		
  for (int i = 0; i < numData; i++)
    data[i] = d(i);
}

FedeasConcr3Material::FedeasConcr3Material(void):
FedeasMaterial(0, MAT_TAG_FedeasConcrete3, 2, 10)
{
	// Does nothing
}

FedeasConcr3Material::~FedeasConcr3Material(void)
{
	// Does nothing
}

UniaxialMaterial*
FedeasConcr3Material::getCopy(void)
{
  Vector d(data, numData);

  FedeasConcr3Material *theCopy = new FedeasConcr3Material(this->getTag(), d);
  
  // Copy history variables
  for (int i = 0; i < 2*numHstv; i++)
    theCopy->hstv[i] = hstv[i];
  
  theCopy->epsilonP = epsilonP;
  theCopy->sigmaP   = sigmaP;
  theCopy->tangentP = tangentP;
  
  return theCopy;
}

double
FedeasConcr3Material::getInitialTangent(void)
{
	//return 2.0*fc/ec;
	return 2.0*data[0]/data[1];
}
