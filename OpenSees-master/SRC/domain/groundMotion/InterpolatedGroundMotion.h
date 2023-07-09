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
// $Date: 2000-12-12 07:32:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/groundMotion/InterpolatedGroundMotion.h,v $
                                                                        
// File: ~/earthquake/InterpolatedGroundMotion.h
// 
// Written: fmk
// Created: 11/00
// Revision: A
//
// Description: This file contains the class definition for 
// an InterpolatedGroundMotion, which inetrploates between 2
// motions to determine a motion.
//
// What: "@(#) InterpolatedGroundMotion.h, revA"

#ifndef InterpolatedGroundMotion_h
#define InterpolatedGroundMotion_h

#include <GroundMotion.h>
#include <TimeSeries.h>
#include <Vector.h>

class InterpolatedGroundMotion : public GroundMotion
{
  public:
    InterpolatedGroundMotion();
    InterpolatedGroundMotion(GroundMotion **groundMotions, 
			     const Vector &fact,
			     bool destroyMotions,
			     double deltaPeak = 0.01);

    virtual ~InterpolatedGroundMotion();

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
    GroundMotion **theMotions;
    Vector *factors;
    int destroyMotions;
    
    Vector data;
    double deltaPeak;  // increment for determining a peak response
};

#endif
