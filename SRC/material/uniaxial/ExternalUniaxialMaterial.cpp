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

// $Revision: 1.0 $
// $Date: 2019-01-26 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/UniaxialMaterial.h,v $

// Written: M. Salehi
// Created: 26-1 2019
//


#include <Vector.h>
#include <string.h>

#include <ExternalUniaxialMaterial.h>
#include <Channel.h>
#include <math.h>
#include <float.h>


ExternalUniaxialMaterial::ExternalUniaxialMaterial(int tag) :UniaxialMaterial(tag, MAT_TAG_ExternalUniaxialMaterial)
{

}

void
ExternalUniaxialMaterial::SetLinks(UMSetTrialStrain _setTrialStrain,
	UMGetStress _getStress,
	UMGetTangent _getTangent,
	UMGetInitialTangent _getInitialTangent,
	UMGetDampTangent _getDampTangent,
	UMGetStrain _getStrain,
	UMGetStrainRate _getStrainRate,
	UMGetRho _getRho,
	UMCommitState _commitState,
	UMRevertToLastCommit _revertToLastCommit,
	UMRevertToStart _revertToStart,
	UMGetCopy _getCopy,
	UMPrint _print) {

	_SetTrialStrain = _setTrialStrain;
	_GetStress = _getStress;
	_GetTangent = _getTangent;
	_GetInitialTangent = _getInitialTangent;
	_GetStrain = _getStrain;
	_CommitState = _commitState;
	_RevertToLastCommit = _revertToLastCommit;
	_RevertToStart = _revertToStart;
	_GetCopy = _getCopy;
	_Print = _print;
	_GetDampTangent = _getDampTangent;
	_GetStrainRate = _getStrainRate;
	_GetRho = _getRho;

}

ExternalUniaxialMaterial::~ExternalUniaxialMaterial()
{
	
}


//int
//ExternalUniaxialMaterial::setTrialStrain(double strain, double strainRate)
//{
//  
//  return 0;
//}

//double
//ExternalUniaxialMaterial::getStress(void)
//{
//  return 0;
//}

//double
//ExternalUniaxialMaterial::getTangent(void)
//{
//  return 0;
//}

//double
//ExternalUniaxialMaterial::getInitialTangent(void)
//{
//  return 0;
//}

//double
//ExternalUniaxialMaterial::getStrain(void)
//{
//  return 0;
//}

//int
//ExternalUniaxialMaterial::commitState(void)
//{
//  return 0;
//}

//int
//ExternalUniaxialMaterial::revertToLastCommit(void)
//{
//  return 0;
//}
//
//int
//ExternalUniaxialMaterial::revertToStart(void)
//{
//  return 0;
//}

//UniaxialMaterial *
//ExternalUniaxialMaterial::getCopy(void)
//{
//  return 0;
//}

int
ExternalUniaxialMaterial::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	return res;
}

int
ExternalUniaxialMaterial::recvSelf(int cTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;
	return res;
}

//void
//ExternalUniaxialMaterial::Print(OPS_Stream &s, int flag)
//{
//	return _PrintSelf(flag);
//}


