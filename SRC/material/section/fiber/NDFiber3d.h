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
                                                                        
// $Revision$
// $Date$
// $Source$


// File: ~/fiber/NDFiber3d.h
//
// Written: MHS
// Created: 2012
// Revision: 
//
// Description: This file contains the class definition for 
// NDFiber3d.h. NDFiber3d provides the abstraction of a
// uniaxial fiber whose position  is defined with only one coordinate.
// The NDFiber3d is subjected to a stress state with 
// only one nonzero axial stress and corresponding axial strain.
//
// What: "@(#) NDFiber3d.h, revA"

#ifndef NDFiber3d_h
#define NDFiber3d_h

#include <Fiber.h>
#include <Vector.h>
#include <Matrix.h>

class NDMaterial;
class Parameter;

class NDFiber3d : public Fiber
{
  public:
    NDFiber3d ();   
    NDFiber3d (int tag, NDMaterial &theMat, double Area, double y, double z);
    ~NDFiber3d();

    
    int   setTrialFiberStrain(const Vector &vs);
    Vector &getFiberStressResultants (void);
    Matrix &getFiberTangentStiffContr (void);

    int   commitState(void);
    int   revertToLastCommit(void);    
    int   revertToStart(void);
    
    Fiber *getCopy(void);
    int getOrder(void);
    const ID &getType(void);

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);
    
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &info);
	
    void getFiberLocation(double &y, double &z);
    NDMaterial *getNDMaterial(void) {return theMaterial;}
    double getArea(void) {return area;};

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int parameterID);

    const Vector &getFiberSensitivity(int gradNumber, bool cond);
    int commitSensitivity(const Vector &dedh, int gradNumber,
			  int numGrads);

  protected:
    
  private:
    NDMaterial *theMaterial;   // pointer to a material
    double area;                          // area of the fiber 
    double y;		// fiber location
    double z;

    static Matrix ks;       // static class wide matrix object for returns
    static Vector fs;	    // static class wide vector object for returns

    static ID code;
};


#endif






