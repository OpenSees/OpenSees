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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/oldUniaxialFiber3d.h,v $
                                                                        
                                                                        
// File: ~/fiber/UniaxialFiber3d.h
//
// Written: Remo Magalhaes de Souza
// Created: 10/98
// Revision: 
//
// Description: This file contains the class definition for 
// UniaxialFiber3d.h. UniaxialFiber3d provides the abstraction of a
// uniaxial fiber that forms a fiber section for 3d frame elements (the 
// fiber position inside the section is defined by two coordinates)
// The UniaxialFiber3d is subjected to a stress state with 
// only one nonzero axial stress and corresponding axial strain.
//
// What: "@(#) UniaxialFiber3d.h, revA"

#ifndef UniaxialFiber3d_h
#define UniaxialFiber3d_h

#include <UniaxialFiber.h>
#include <Matrix.h>
#include <Vector.h>

class UniaxialMaterial;

class UniaxialFiber3d: public UniaxialFiber
{
  public:
    UniaxialFiber3d (int tag, UniaxialMaterial &theMat, double Area, 
                     const Vector &position);
        UniaxialFiber3d ();
 
    ~UniaxialFiber3d();

    int   setTrialFiberStrain(const Vector &vs);
    Vector &getFiberStressResultants (void);
    Matrix &getFiberTangentStiffContr (void); 
    int   commitState(void);
    int   revertToLastCommit(void);    
    int   revertToStart(void);
    
    UniaxialFiber3d *getCopy(void);

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(ostream &s, int flag =0);

    Vector getPosition (void) const;
    double getArea     (void) const;
    double getStress   (void) const;
    double getStrain   (void) const;
         
  protected:
    
  private:
    UniaxialMaterial *theMaterial;  // pointer to a material
    double area;                          // area of the fiber 
    Matrix as;                            // matrix that transforms
                           // section deformations into fiber strain
    Matrix ks;
    Vector fs;
    
};


#endif

