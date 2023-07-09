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
                                                                        
// $Revision: 1.6 $
// $Date: 2008-03-05 20:21:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/groundMotion/GroundMotion.h,v $
                                                                        
                                                                        
// File: ~/earthquake/GroundMotion.h
// 
// Written: fmk 
// Created: 05/98
// Revision: A
//
// Description: This file contains the class definition for GroundMotion.
// GroundMotion is an abstract base class, i.e. no objects of it's
// type can be created. A GroundMotion describes the motion of the ground
// in a dynamic analysis.
//
// What: "@(#) GroundMotion.h, revA"

#ifndef GroundMotion_h
#define GroundMotion_h

#include <MovableObject.h>
#include <TimeSeries.h>
#include <TimeSeriesIntegrator.h>
#include <Vector.h>

class GroundMotion : public MovableObject
{
  public:
    GroundMotion(TimeSeries *dispSeries, 
		 TimeSeries *velSeries, 
		 TimeSeries *accelSeries,
		 TimeSeriesIntegrator *theIntegrator = 0,
		 double dTintegration = 0.01, double fact = 1.0);
    
    GroundMotion(int classTag);
    virtual ~GroundMotion();

    virtual double getDuration(void);

    virtual double getPeakAccel();
    virtual double getPeakVel();
    virtual double getPeakDisp();

    virtual double getAccel(double time);
    virtual double getVel(double time);
    virtual double getDisp(double time);
    virtual const  Vector &getDispVelAccel(double time);
    
    void setIntegrator(TimeSeriesIntegrator *integrator);
    TimeSeries *integrate(TimeSeries *theSeries, double delta = 0.01); 

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    const TimeSeries *getAccelSeries(void) const {return theAccelSeries;}
	
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    virtual double getAccelSensitivity(double time);
    virtual int setParameter(const char **argv, int argc, Parameter &param);
    //    virtual int updateParameter(int parameterID, Information &info);
    // virtual int activateParameter(int parameterID);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:

  private:
    TimeSeries *theAccelSeries;  // Ground acceleration
    TimeSeries *theVelSeries;	 // Ground velocity
    TimeSeries *theDispSeries;	 // Ground displacement
    TimeSeriesIntegrator *theIntegrator;
    
    Vector data;
    double delta;
    double fact;
};

#endif
