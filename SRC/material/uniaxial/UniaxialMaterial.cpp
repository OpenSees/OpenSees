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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/UniaxialMaterial.cpp,v $
                                                                        
                                                                        
// File: ~/material/UniaxialMaterial.C
//
// Written: fmk 
// Created: 05/98
// Revision: A
//
// Description: This file contains the class implementation for 
// UniaxialMaterial.
//
// What: "@(#) UniaxialMaterial.C, revA"

#include <UniaxialMaterial.h>
#include <string.h>
#include <Information.h>

UniaxialMaterial::UniaxialMaterial(int tag, int clasTag)
:Material(tag,clasTag)
{

}

UniaxialMaterial::~UniaxialMaterial()
{
  // does nothing
}

// default operation for strain rate is zero
double
UniaxialMaterial::getStrainRate(void)
{
    return 0.0;
}

// default operation for damping tangent is zero
double
UniaxialMaterial::getDampTangent(void)
{
    return 0.0;
}

// default operation for secant stiffness is the tangent
double
UniaxialMaterial::getSecant (void)
{
    return this->getTangent();
}

UniaxialMaterial*
UniaxialMaterial::getCopy(SectionForceDeformation *s)
{
	return this->getCopy();
}

int 
UniaxialMaterial::setResponse(char **argv, int argc, Information &matInfo)
{
    // stress
    if (strcmp(argv[0],"stress") ==0) {
	matInfo.theType = DoubleType;
	return 1;
    } 

    // tangent
    else if (strcmp(argv[0],"tangent") ==0) {
	matInfo.theType = DoubleType;
	return 2;
    }

    // otherwise unknown
    else
	return -1;
}

int 
UniaxialMaterial::getResponse(int responseID, Information &matInfo)
{
    // each subclass must implements it's own stuff    
  switch (responseID) {
    case -1:
      return -1;
      
    case 1:
      matInfo.theDouble = this->getStress();    
      return 0;
      
    case 2:
      matInfo.theDouble = this->getTangent();    
      return 0;      
      
    default:      
      return -1;
  }
    
}
