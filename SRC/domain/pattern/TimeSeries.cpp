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
// $Date: 2003-03-04 00:48:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TimeSeries.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for TimeSeries.
// TimeSeries is an abstract class. A TimeSeries object is used to
// determine the load factor to be applied to the loads in a pettern.
// to the model. 
//
// What: "@(#) TimeSeries.C, revA"

#include <TimeSeries.h>

TimeSeries::TimeSeries(int classTag)
  :MovableObject(classTag)
{

}

TimeSeries::~TimeSeries()
{

}

// AddingSensitivity:BEGIN //////////////////////////////////////////
int
TimeSeries::setParameter(const char **argv, int argc, Information &info)
{
	return 0;
}

int
TimeSeries::updateParameter(int parameterID, Information &info)
{
	return 0;
}

int
TimeSeries::activateParameter(int parameterID)
{
	return 0;
}

double
TimeSeries::getFactorSensitivity(double pseudoTime)
{
	return 0.0;
}

// AddingSensitivity:END ////////////////////////////////////////////
