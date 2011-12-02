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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-25 23:33:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/Material.cpp,v $
                                                                        
                                                                        
// File: ~/material/Material.C
//
// Written: fmk 
// Created: 05/98
// Revision: A
//
// Description: This file contains the class implementation for MaterialModel.
//
// What: "@(#) MaterialModel.C, revA"

#include <Material.h>

Material::Material(int tag, int clasTag)
:TaggedObject(tag), MovableObject(clasTag)
{

}


Material::~Material()
{
  // does nothing


}

int
Material::setVariable(const char *argv)
{
	return -1;
}

int
Material::getVariable(int variableID, double &info)
{
	return -1;
}

int
Material::setParameter(const char **argv, int argc, Information &eleInformation)
{
    return -1;
}

int
Material::updateParameter(int responseID, Information &eleInformation)
{
    return -1;
}

Response*
Material::setResponse(const char **argv, int argc, Information &info)
{
	return 0;
}

int 
Material::getResponse(int responseID, Information &info)
{
	return -1;
}


