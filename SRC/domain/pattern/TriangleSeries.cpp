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

// $Revision: 1.1 $
// $Date: 2005-12-15 00:35:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TriangleSeries.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/04
// Revision: A
//
// Purpose: This file contains the class implementation of TriangleSeries.
//
// What: "@(#) TriangleSeries.cpp, revA"


#include <TriangleSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>

#include <math.h>


TriangleSeries::TriangleSeries(double startTime, double finishTime,
    double T, double phi, double theFactor)
    : TimeSeries(TSERIES_TAG_TriangleSeries),
    tStart(startTime),tFinish(finishTime),
    period(T),shift(phi),cFactor(theFactor)
{
	if (period == 0.0)  {
		opserr << "TriangleSeries::TriangleSeries -- input period is zero, setting period to 1\n";
		period = 1;
	}
}


TriangleSeries::TriangleSeries()
    : TimeSeries(TSERIES_TAG_TriangleSeries),
    tStart(0.0),tFinish(0.0),
    period(1.0),shift(0.0),cFactor(1.0)
{
	// does nothing
}


TriangleSeries::~TriangleSeries()
{
	// does nothing
}


double TriangleSeries::getFactor(double pseudoTime)
{	
	if (tStart <= pseudoTime && pseudoTime <= tFinish)  {
		double slope = cFactor/(period/4);
		double k = (pseudoTime+shift)/period - floor((pseudoTime+shift)/period);
		if (k < 0.25)
			return slope*k*period;
		else if (k < 0.75)
			return cFactor-slope*(k-0.25)*period;
		else if (k < 1.00)
			return -cFactor + slope*(k-0.75)*period;
		else
			return 0;
	}
	else
		return 0;
}


int TriangleSeries::sendSelf(int commitTag, Channel &theChannel)
{
	int dbTag = this->getDbTag();
	Vector data(5);
	data(0) = cFactor;
	data(1) = tStart;	
	data(2) = tFinish;
	data(3) = period;
	data(4) = shift;
	int result = theChannel.sendVector(dbTag,commitTag, data);
	if (result < 0)  {
		opserr << "TriangleSeries::sendSelf() - channel failed to send data\n";
		return result;
	}
	return 0;
}


int TriangleSeries::recvSelf(int commitTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
	int dbTag = this->getDbTag();
	Vector data(5);
	int result = theChannel.recvVector(dbTag,commitTag, data);
	if (result < 0)  {
		opserr << "TriangleSeries::sendSelf() - channel failed to receive data\n";
		cFactor = 1.0;
		tStart  = 0.0;
		tFinish = 0.0;
		period  = 1.0;
		shift   = 0.0;
		return result;
	}
	cFactor = data(0);
	tStart  = data(1);
	tFinish = data(2);
	period  = data(3);
	shift   = data(4);
	
	return 0;    
}


void TriangleSeries::Print(OPS_Stream &s, int flag)
{
	s << "Triangle Series" << endln;
	s << "\tFactor: " << cFactor << endln;
	s << "\ttStart: " << tStart << endln;
	s << "\ttFinish: " << tFinish << endln;
	s << "\tPeriod: " << period << endln;
	s << "\tPhase Shift: " << shift << endln;
}
