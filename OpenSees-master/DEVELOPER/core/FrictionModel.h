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

// $Revision: 4952 $
// $Date: 2012-08-08 22:56:05 -0700 (Wed, 08 Aug 2012) $
// $URL: svn://opensees.berkeley.edu/usr/local/svn/OpenSees/trunk/SRC/element/frictionBearing/frictionModel/FrictionModel.h $

#ifndef FrictionModel_h
#define FrictionModel_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class definition for FrictionModel.
// FrictionModel is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 

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
    virtual double getNormalForce();
    virtual double getVelocity();
    virtual double getFrictionForce() = 0;
    virtual double getFrictionCoeff() = 0;
    virtual double getDFFrcDNFrc() = 0;
    virtual double getDFFrcDVel() = 0;
    
    virtual int commitState() = 0;
    virtual int revertToLastCommit() = 0;
    virtual int revertToStart() = 0;
    
    virtual FrictionModel *getCopy() = 0;
    
    virtual Response *setResponse(const char **argv, int argc,
        OPS_Stream &theOutputStream);
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

extern bool OPS_addFrictionModel(FrictionModel *newComponent);
extern FrictionModel *OPS_getFrictionModel(int tag);
extern void OPS_clearAllFrictionModel();

#endif
