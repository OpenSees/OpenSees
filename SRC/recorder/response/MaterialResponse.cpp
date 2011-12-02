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
// $Date: 2000-12-18 11:35:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/MaterialResponse.cpp,v $
                                                                        
// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the MaterialResponse class implementation

#include <MaterialResponse.h>
#include <Material.h>

MaterialResponse::MaterialResponse(Material *mat, int id):
Response(), theMaterial(mat), responseID(id), matInfo()
{

}

MaterialResponse::MaterialResponse(Material *mat, int id, int val):
Response(), theMaterial(mat), responseID(id), matInfo(val)
{

}

MaterialResponse::MaterialResponse(Material *mat, int id, double val):
Response(), theMaterial(mat), responseID(id), matInfo(val)
{

}

MaterialResponse::MaterialResponse(Material *mat, int id, const ID &val):
Response(), theMaterial(mat), responseID(id), matInfo(val)
{

}

MaterialResponse::MaterialResponse(Material *mat, int id, const Vector &val):
Response(), theMaterial(mat), responseID(id), matInfo(val)
{

}

MaterialResponse::MaterialResponse(Material *mat, int id, const Matrix &val):
Response(), theMaterial(mat), responseID(id), matInfo(val)
{

}

MaterialResponse::MaterialResponse(Material *mat, int id, const Tensor &val):
Response(), theMaterial(mat), responseID(id), matInfo(val)
{

}

MaterialResponse::~MaterialResponse()
{

}

int
MaterialResponse::getResponse(void)
{
	return theMaterial->getResponse(responseID, matInfo);
}

void
MaterialResponse::Print(ostream &s, int flag)
{
	matInfo.Print(s, flag);
}
