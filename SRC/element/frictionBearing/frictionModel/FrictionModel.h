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
// $Date: 2009-04-17 23:02:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/frictionBearing/frictionModel/FrictionModel.h,v $

#ifndef FrictionModel_h
#define FrictionModel_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class definition for FrictionModel.
// FrictionModel is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// What: "@(#) FrictionModel.h, revA"

#include <DomainComponent.h>
#include <MovableObject.h>
#include <TaggedObject.h>
#include <Vector.h>

class Response;

class FrictionModel : public TaggedObject, public MovableObject
{
public:
    // constructor
    FrictionModel(int tag, int classTag);
    
    // destructor
    virtual ~FrictionModel();
    
    // public methods to set and obtain response
    virtual int setTrial(double normalForce, double velocity = 0.0) = 0;
    virtual double getNormalForce(void);
    virtual double getVelocity(void);
    virtual double getFrictionForce(void) = 0;
    virtual double getFrictionCoeff(void) = 0;
    virtual double getDFFrcDNFrc(void) = 0;
    
    virtual int commitState(void) = 0;
    virtual int revertToLastCommit(void) = 0;
    virtual int revertToStart(void) = 0;
    
    virtual FrictionModel *getCopy(void) = 0;
    
    virtual Response *setResponse(char **argv, int argc, Information &info);
    virtual int getResponse(int responseID, Information &info);
    
    virtual int sendSelf(int commitTag, Channel &theChannel) = 0;  
    virtual int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker) = 0;
    
    virtual void Print(OPS_Stream &s, int flag = 0) = 0;
    
protected:
    double trialN;      // trial normal contact force
    double trialVel;    // trial velocity
    
private:

};

#endif
