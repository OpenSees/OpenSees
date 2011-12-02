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
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/UniaxialFiber2d.h,v $
                                                                        
                                                                        
// File: ~/fiber/UniaxialFiber2d.h
//
// Written: Remo Magalhaes de Souza
// Created: 10/98
// Revision: 
//
// Description: This file contains the class definition for 
// UniaxialFiber2d.h. UniaxialFiber2d provides the abstraction of a
// uniaxial fiber whose position  is defined with only one coordinate.
// The UniaxialFiber2d is subjected to a stress state with 
// only one nonzero axial stress and corresponding axial strain.
//
// What: "@(#) UniaxialFiber2d.h, revA"

#ifndef UniaxialFiber2d_h
#define UniaxialFiber2d_h

#include <Fiber.h>
#include <Vector.h>
#include <Matrix.h>

class UniaxialMaterial;

class UniaxialFiber2d : public Fiber
{
  public:
    UniaxialFiber2d ();   
    UniaxialFiber2d (int tag, UniaxialMaterial &theMat, double Area, double position);
    ~UniaxialFiber2d();

    
    int   setTrialFiberStrain(const Vector &vs);
    Vector &getFiberStressResultants (void);
    Matrix &getFiberTangentStiffContr (void);
    Matrix &getFiberSecantStiffContr (void);

    int   commitState(void);
    int   revertToLastCommit(void);    
    int   revertToStart(void);
    
    Fiber *getCopy(void);
	int getOrder(void);
	const ID &getType(void);

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(ostream &s, int flag =0);
    
  protected:
    
  private:
    UniaxialMaterial *theMaterial;   // pointer to a material
    double area;                          // area of the fiber 
	double y;		// fiber location

    static Matrix ks;       // static class wide matrix object for returns
    static Vector fs;	    // static class wide vector object for returns

	static ID code;
};


#endif






