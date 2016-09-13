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
// $Date: 2010-02-04 00:34:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/ConstantSeries.h,v $

#ifndef ConstantSeries_h
#define ConstantSeries_h

// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for ConstantSeries.
// ConstantSeries is a concrete class. A ConstantSeries object provides
// a linear time series. the factor is given by the pseudoTime and 
// a constant factor provided in the constructor.
//
// What: "@(#) ConstantSeries.h, revA"

#include <TimeSeries.h>

class ConstantSeries : public TimeSeries
{
  public:
    ConstantSeries(int tag=0, double cFactor = 1.0);

    ~ConstantSeries();

    TimeSeries *getCopy(void);

    // method to get load factor
    double getFactor(double pseudoTime) {return cFactor;}

    // None of the following functions should be invoked on this type
    // of object
    double getDuration () {return 0.0;} // dummy function
    double getPeakFactor () {return cFactor;} // dummy function
    double getTimeIncr (double pseudoTime) {return 1.0;} // dummy function
    
    // methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    double getFactorSensitivity(double pseudoTime);
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int parameterID);
    // AddingSensitivity:BEGIN //////////////////////////////////////////

  protected:
	
  private:
    double cFactor;  // factor = cFactor
    int parameterID;
};

#endif
