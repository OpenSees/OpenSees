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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/oldUniaxialFiber3d.cpp,v $
                                                                        
                                                                        
// File: ~/section/UniaxialFiber3d.C
//
// Written: Remo Magalhaes de Souza
// Created: 10/98
// Revision: 
//
// Description: This file contains the implementation for the
// UniaxialFiber3d class. UniaxialFiber3d provides the abstraction of a
// uniaxial fiber that forms a fiber section for 3d frame elements.
// The UniaxialFiber3d is subjected to a stress state with
// only one nonzero axial stress and corresponding axial strain.
//
// What: "@(#) UniaxialFiber3d.C, revA"

#include <stdlib.h>

#include <UniaxialMaterial.h>
#include <UniaxialFiber3d.h>
#include <Vector.h>

// constructor:
UniaxialFiber3d::UniaxialFiber3d(int tag, 
                                 UniaxialMaterial &theMat,
                                 double Area, const Vector &position):
                                 UniaxialFiber(tag, FIB_TAG_UniaxialFiber3d),
                                 theMaterial(0), area(Area), as(1,3), fs(3), ks(3,3)
{
   theMaterial = theMat.getCopy();  // get a copy of the MaterialModel
   as(0,0) = 1.0;
   as(0,1) = -position(0);
   as(0,2) =  position(1);
}

UniaxialFiber3d::UniaxialFiber3d()
:	
UniaxialFiber(0, FIB_TAG_UniaxialFiber3d),
theMaterial(0), area(0), as(1,3), fs(3), ks(3,3)
{
   as(0,0) = 0;
   as(0,1) = 0;
   as(0,2) = 0;
}

// destructor:
UniaxialFiber3d::~UniaxialFiber3d ()
{
   if (theMaterial != 0)
      delete theMaterial;
}


int   
UniaxialFiber3d::setTrialFiberStrain(const Vector &vs)
{
   Vector strain(1);
   strain = as * vs;
  
   theMaterial->setTrialStrain (strain(0));

   return 0;
}



// get fiber stress resultants 
Vector &
UniaxialFiber3d::getFiberStressResultants (void)
{
    double df;              // axial force in fiber
    double stress;          // stress in fiber
//    Vector fs(3);           // stress resultants 
   
    stress = theMaterial->getStress();  // get the stress from the 
                                        // material 
    df = stress * area;

//    fs = as^ df;
   for (int i = 0; i < 3; i++)
     fs(i) = as(0,i) * df;

    return fs;
}



// get contribution of fiber to section tangent stiffness
Matrix &
UniaxialFiber3d::getFiberTangentStiffContr(void) 
{
    double Et;              // material tangent
//    Matrix ks(3,3);         // fiber tangent stiff 
   
    Et = theMaterial->getTangent();  // get the material tangent

    ks = (as^as) * area * Et;      // as^as = asT  * as

    return ks;
}



UniaxialFiber3d *
UniaxialFiber3d::getCopy (void)
{
   // make a copy of the fiber 
   Vector position(2);

   position(0) = -as(0,1);
   position(1) =  as(0,2);

   UniaxialFiber3d *theCopy = new UniaxialFiber3d (this->getTag(), 
                                                   *theMaterial, area, 
                                                   position);
   return theCopy;
}  


int   
UniaxialFiber3d::commitState(void)
{
   return theMaterial->commitState();
}


int   
UniaxialFiber3d::revertToLastCommit(void)
{
   return theMaterial->revertToLastCommit();
}

int   
UniaxialFiber3d::revertToStart(void)
{
   return theMaterial->revertToStart();
}


int   
UniaxialFiber3d::sendSelf(int cTag, Channel &theChannel)
{
   // TO DO
   cerr << "FATAL - UniaxialFiber3d::sendSelf () - not yet implemented";
   exit (-1);
   return 0;
}


int   
UniaxialFiber3d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
   // TO DO
   cerr << "FATAL - UniaxialFiber3d::recvSelf () - not yet implemented";
   exit (-1);
   return 0;
}


Vector  
UniaxialFiber3d::getPosition(void) const
{
   Vector position(2);
   position(0) = -as(0,1);
   position(1) =  as(0,2);

   return position;
}

double 
UniaxialFiber3d::getArea(void) const
{
   return area;
}


double 
UniaxialFiber3d::getStress(void) const
{
   return theMaterial->getStress();
}


double 
UniaxialFiber3d::getStrain(void) const
{
   return theMaterial->getStrain();
}


void UniaxialFiber3d::Print(ostream &s, int flag)
{
   if (flag == 1)
   {
      Vector position = this->getPosition();

      double strain = theMaterial->getStrain();
      double stress = theMaterial->getStress();

      s << "\n" << position(0) << " " << position(1) << " "  << area  << " " << stress << " " << strain;
 
   }
   else
   {

//    s << "\nFiber: " << this->getTag() << " Type: UniaxialFiber3d";
//    s << "\tMaterial: " << *theMaterial;
//    s << "\tMatrix as:\n" << as; 
//    s << "\tArea: " << area; 

    char st[200];
    sprintf (st, "Material=  %3d  y=%12.4e  z=%12.4e  A=%12.4e\n",theMaterial->getTag(),
             -as(0,1),as(0,2),area);   

//    s << " Material: " << theMaterial->getTag();
//    s << "  y: " << -as(0,1) << " z: " <<  as(0,2);
//    s << "  Area: " << area << endl; 
      s << st;
   }
}
   
   
ostream &operator<<(ostream &s, UniaxialFiber3d &E)  
{
    E.Print(s); 
    return s; 
}

