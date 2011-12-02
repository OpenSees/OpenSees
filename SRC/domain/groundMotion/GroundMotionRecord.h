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
// $Date: 2003-02-25 23:32:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/groundMotion/GroundMotionRecord.h,v $
                                                                        
                                                                        
// File: ~/earthquake/GroundMotionRecord.h
// 
// Written: MHS
// Created: 09/99
// Revision: A
//
// Description: This file contains the class definition for 
// a GroundMotionRecord, which creates a GroundMotion from
// user supplied data contained in an input file.
//
// What: "@(#) GroundMotionRecord.h, revA"

#ifndef GroundMotionRecord_h
#define GroundMotionRecord_h

#include <GroundMotion.h>
#include <TimeSeries.h>
#include <Vector.h>

class GroundMotionRecord : public GroundMotion
{
  public:
    GroundMotionRecord();
    GroundMotionRecord(const char *fileNameAccel, double timeStep,
		       double fact = 1.0, double dTintegration = 0.01);
    GroundMotionRecord(const char *fileNameAccel, const char *fileNameTime,
		       double fact = 1.0, double dTintegration = 0.01);
    // ADD CONSTRUCTORS FOR OTHER INPUT FORMATS

    virtual ~GroundMotionRecord();

    virtual double getDuration();

    virtual double getPeakAccel();
    virtual double getPeakVel();
    virtual double getPeakDisp();

    virtual double getAccel(double time);    
    virtual double getVel(double time);    
    virtual double getDisp(double time);        
    virtual const  Vector &getDispVelAccel(double time);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
  protected:

  private:
    TimeSeries *theAccelTimeSeries;  // Ground acceleration
    TimeSeries *theVelTimeSeries;	 // Ground velocity
    TimeSeries *theDispTimeSeries;	 // Ground displacement
    
    Vector data;
    double delta;
};

#endif
