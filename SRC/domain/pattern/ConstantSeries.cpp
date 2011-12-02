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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/ConstantSeries.cpp,v $

// File: ~/domain/pattern/ConstantSeries.C
//
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for ConstantSeries.
// ConstantSeries is a concrete class. A ConstantSeries object provides
// a linear time series. the factor is given by the pseudoTime and 
// a constant factor provided in the constructor.
//
// What: "@(#) ConstantSeries.C, revA"


#include <ConstantSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>

ConstantSeries::ConstantSeries(double theFactor)
  :TimeSeries(TSERIES_TAG_ConstantSeries),
   cFactor(theFactor)
{
  // does nothing
}


ConstantSeries::~ConstantSeries()
{
  // does nothing
}

int
ConstantSeries::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  Vector data(1);
  data(0) = cFactor;
  int result = theChannel.sendVector(dbTag,commitTag, data);
  if (result < 0) {
    cerr << "ConstantSeries::sendSelf() - channel failed to send data\n";
    return result;
  }
  return 0;
}


int 
ConstantSeries::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  Vector data(1);
  int result = theChannel.recvVector(dbTag,commitTag, data);
  if (result < 0) {
    cerr << "ConstantSeries::sendSelf() - channel failed to receive data\n";
    cFactor = 1.0;
    return result;
  }
  cFactor = data(0);

  return 0;    
}


void
ConstantSeries::Print(ostream &s, int flag)
{
    s << "Constant Series: factor: " << cFactor << "\n";

}
