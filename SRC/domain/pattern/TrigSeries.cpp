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
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TrigSeries.cpp,v $
                                                                        
                                                                        
// File: ~/domain/pattern/TrigSeries.C
//
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class implementation of TrigSeries.
//
// What: "@(#) TrigSeries.C, revA"


#include <TrigSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>

#include <math.h>

TrigSeries::TrigSeries(double startTime, double finishTime,
					   double T, double phi, double theFactor)
  :TimeSeries(TSERIES_TAG_TrigSeries),
   tStart(startTime),tFinish(finishTime),
   period(T),shift(phi),cFactor(theFactor)
{
	if (period == 0.0) {
		g3ErrorHandler->warning("%s -- input period is zero, setting period to PI",
		"TrigSeries::TrigSeries");
		period = 2*asin(1.0);
	}
}

TrigSeries::TrigSeries()
  :TimeSeries(TSERIES_TAG_TrigSeries),
   tStart(0.0),tFinish(0.0),period(1.0),shift(0.0),cFactor(1.0)
{
  // does nothing
}


TrigSeries::~TrigSeries()
{
  // does nothing
}

double
TrigSeries::getFactor(double pseudoTime)
{
	static double twopi = 4*asin(1.0);

	if (pseudoTime >= tStart && pseudoTime <= tFinish)
		return cFactor*sin(twopi*(pseudoTime-tStart)/period + shift);
	else
		return 0.0;
}

int
TrigSeries::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  Vector data(5);
  data(0) = cFactor;
  data(1) = tStart;	
  data(2) = tFinish;
  data(3) = period;
  data(4) = shift;
  int result = theChannel.sendVector(dbTag,commitTag, data);
  if (result < 0) {
    cerr << "TrigSeries::sendSelf() - channel failed to send data\n";
    return result;
  }
  return 0;
}


int 
TrigSeries::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  Vector data(5);
  int result = theChannel.recvVector(dbTag,commitTag, data);
  if (result < 0) {
    cerr << "TrigSeries::sendSelf() - channel failed to receive data\n";
    cFactor = 1.0;
    tStart= 0.0;
    tFinish = 0.0;
	period = 1.0;
	shift = 0.0;
    return result;
  }
  cFactor = data(0);
  tStart = data(1);
  tFinish = data(2);
  period = data(3);
  shift = data(4);

  return 0;    
}

void
TrigSeries::Print(ostream &s, int flag)
{
    s << "Trig Series" << endl;
	s << "\tFactor: " << cFactor << endl;
	s << "\ttStart: " << tStart << endl;
	s << "\ttFinish: " << tFinish << endl;
	s << "\tPeriod: " << period << endl;
	s << "\tPhase Shift: " << shift << endl;
}
