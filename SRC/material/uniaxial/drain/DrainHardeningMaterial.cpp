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
                                                                        
// $Revision: 1.2 $
// $Date: 2001-10-01 17:08:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/drain/DrainHardeningMaterial.cpp,v $
                                                                      
// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// DrainHardeningMaterial.

#include <DrainHardeningMaterial.h>
#include <Vector.h>

DrainHardeningMaterial::DrainHardeningMaterial(int tag,
	double E, double sigY, double Hiso, double Hkin, double b):
// 3 history variables and 4 material parameters
DrainMaterial(tag, MAT_TAG_DrainHardening, 3, 4, b)
{
	data[0] = E;
	data[1] = sigY;
	data[2] = Hiso;
	data[3] = Hkin;
}

DrainHardeningMaterial::DrainHardeningMaterial(void):
DrainMaterial(0, MAT_TAG_DrainHardening, 3, 4)
{
	// Does nothing
}

DrainHardeningMaterial::~DrainHardeningMaterial(void)
{
	// Does nothing
}
