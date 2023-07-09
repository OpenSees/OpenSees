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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/drain/DrainClough1Material.cpp,v $

// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// DrainClough1Material.

#include <DrainClough1Material.h>
#include <Vector.h>

DrainClough1Material::DrainClough1Material(int tag,
	double E, double fyp, double fyn, double alpha,
	double ecaps, double ecapk, double ecapa, double ecapd,
	double cs, double ck, double ca, double cd,
	double capSlope, double capDispP, double capDispN, double res,
	double b):
// 15 history variables and 16 material parameters
DrainMaterial(tag, MAT_TAG_DrainClough1, 15, 16, b)
{
	data[0]  = E;
	data[1]  = fyp;
	data[2]  = fyn;
	data[3]  = alpha;
	data[4]  = ecaps;
	data[5]  = ecapk;
	data[6]  = ecapa;
	data[7]  = ecapd;
	data[8]  = cs;
	data[9]  = ck;
	data[10] = ca;
	data[11] = cd;
	data[12] = capSlope;
	data[13] = capDispP;
	data[14] = capDispN;
	data[15] = res;

	// Initialize history variables
	this->revertToStart();
}

DrainClough1Material::DrainClough1Material(int tag, const Vector &input, double b):
// 15 history variables and 16 material parameters
DrainMaterial(tag, MAT_TAG_DrainClough1, 15, 16, b)
{
	for (int i = 0; i < 16; i++)
		data[i] = input(i);

	// Initialize history variables
	this->revertToStart();
}

DrainClough1Material::DrainClough1Material(void):
DrainMaterial(0, MAT_TAG_DrainClough1, 15, 16)
{
	// Does nothing
}

DrainClough1Material::~DrainClough1Material(void)
{
	// Does nothing
}

int
DrainClough1Material::revertToStart(void)
{
	double dyp = data[1]/data[0];	// fyp/E
	double dyn = data[2]/data[0];	// fyn/E

	hstv[0]  = data[0];		// E
	hstv[1]  = data[0];		// E
	hstv[2]  = dyp;
	hstv[3]  = dyn; 
	hstv[4]  = 0.0;
	hstv[5]  = dyp;
	hstv[6]  = dyn;
	hstv[7]  = data[1];		// fyp
	hstv[8]  = data[2];		// fyn
	hstv[9]  = data[13];	// capDispP
	hstv[10] = data[14];	// capDispN
	hstv[11] = 0.0;
	hstv[12] = 0.0;
	hstv[13] = 0.0;
	hstv[14] = data[0];		// E

	// Set trial history variables to committed values
	for (int i = 0; i < 15; i++)
		hstv[i+15] = hstv[i];

	return 0;
}

UniaxialMaterial*
DrainClough1Material::getCopy(void)
{
	Vector input(data, 16);

	DrainClough1Material *theCopy =
		new DrainClough1Material(this->getTag(), input, beto);

	return theCopy;
}
