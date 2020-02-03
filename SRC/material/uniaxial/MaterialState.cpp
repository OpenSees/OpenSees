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
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: MHS
// Created: August 2000
//
// Description: This file contains the interface for MaterialState,
// which is the base for classes which model states of hysteretic
// degradation.

#include <MaterialState.h>

MaterialState::MaterialState(int tag, int classTag)
  :TaggedObject(tag), MovableObject(classTag)
{
  
}

MaterialState::~MaterialState()
{
  
}

int
MaterialState::setVariable(const char *argv)
{
  return -1;
}

int
MaterialState::getVariable(int variableID, double &info)
{
  return -1;
}

int
MaterialState::setParameter(const char **argv, int argc,
			    Information &eleInformation)
{
  return -1;
}

int
MaterialState::updateParameter(int responseID, Information &eleInformation)
{
  return -1;
}
