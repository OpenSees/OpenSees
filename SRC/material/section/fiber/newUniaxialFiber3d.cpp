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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/newUniaxialFiber3d.cpp,v $
                                                                        
                                                                        
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
#include <stdio.h>

#include <UniaxialMaterial.h>
#include <UniaxialFiber3d.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>

Matrix UniaxialFiber3d::ks(3,3); 
Vector UniaxialFiber3d::fs(3); 

// constructor:
UniaxialFiber3d::UniaxialFiber3d()
:UniaxialFiber(0, FIB_TAG_UniaxialFiber3d),
 theMaterial(0), area(0.0)
{
   as[0] = 0.0;
   as[1] = 0.0;
}

UniaxialFiber3d::UniaxialFiber3d(int tag, 
                                 UniaxialMaterial &theMat,
                                 double Area, const Vector &position)
:UniaxialFiber(tag, FIB_TAG_UniaxialFiber3d),
 theMaterial(0), area(Area)
{
   theMaterial = theMat.getCopy();  // get a copy of the MaterialModel
   as[0] = -position(0);
   as[1] =  position(1);
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
  // Explicitly carry out the multiplication rather than calling
  // the overloaded matrix multiply operator
  double strain = vs(0) + as[0]*vs(1) + as[1]*vs(2);

  if (theMaterial != 0)
      theMaterial->setTrialStrain(strain);
  else {
      g3ErrorHandler->fatal("UniaxialFiber3d::setTrialFiberStrain() - %s\n",
			    "recvSelf() has not been invoked");
      return -1; // in case fatal does not exit
  }

   return 0;
}



// get fiber stress resultants 
Vector &
UniaxialFiber3d::getFiberStressResultants (void) const
{
    double df;              // axial force in fiber
    double stress;          // stress in fiber
   
    stress = theMaterial->getStress();  // get the stress from the 
                                        // material 
    df = stress * area;

    // fs = as^ df;
    fs(0) = df;
    fs(1) = as[0]*df;
    fs(2) = as[1]*df;

    return fs;
}



// get contribution of fiber to section tangent stiffness
Matrix &
UniaxialFiber3d::getFiberTangentStiffContr(void) const
{
    double Et;    // material tangent
   
    Et = theMaterial->getTangent();  // get the material tangent

    // ks = (as^as) * area * Et;      // as^as = asT  * as
    //  which is written as:
    double value = area*Et;
    double as1 = as[0];
    double as2 = as[1];
    double vas1 = as1*value;
    double vas2 = as2*value;
    double vas1as2 = vas1*as2;

    ks(0,0) = value;
    ks(0,1) = vas1;
    ks(0,2) = vas2;
    
    ks(1,0) = ks(0,1);
    ks(1,1) = vas1*as1;
    ks(1,2) = vas1as2;
    
    ks(2,0) = vas2;
    ks(2,1) = vas1as2;
    ks(2,2) = vas2*as2;

    return ks;
}


Matrix &
UniaxialFiber3d::getFiberSecantStiffContr (void) const
{
    // Material secant stiffness
    double S;

    // Get the material secant stiffness
    S = theMaterial->getSecant ();

    // Use the section kinematic matrix to get the 
    // fiber secant stiffness matrix

    // ks = (as^as) * area * S;      // as^as = asT  * as
    //  which is written as:
    double value = area*S;
    double as1 = as[0];
    double as2 = as[1];
    double vas1 = as1*value;
    double vas2 = as2*value;
    ks(0,0) = value;
    ks(0,1) = vas1;
    ks(0,2) = vas2;
    
    ks(1,0) = vas1;
    ks(1,1) = vas1*as1;
    ks(1,2) = vas1*as2;
    
    ks(2,0) = vas2;
    ks(2,1) = vas2*as1;
    ks(2,2) = vas2*as2;
 
    return ks;
}


UniaxialFiber3d *
UniaxialFiber3d::getCopy (void)
{
   // make a copy of the fiber 
   Vector position(2);

   position(0) = -as[0];
   position(1) =  as[1];

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
UniaxialFiber3d::sendSelf(int commitTag, Channel &theChannel)
{
    // 
    // store tag and material info in an ID and send it
    //

    static ID idData(3);
    int dbTag = this->getDbTag();
    idData(0) = this->getTag();
    idData(1) = theMaterial->getClassTag();
    int matDbTag = theMaterial->getDbTag();
    if (matDbTag == 0) {
	matDbTag = theChannel.getDbTag();
	if (matDbTag != 0)
	    theMaterial->setDbTag(matDbTag);
    }
    idData(2) = matDbTag;
    
    if (theChannel.sendID(dbTag, commitTag, idData) < 0)  {
	g3ErrorHandler->warning("UniaxialFiber3d::sendSelf() - %s\n",
			      "failed to send ID data");
	return -1;
    }    
    
    // 
    // store area and position data in a vector and send it
    //
    
    static Vector dData(3);
    dData(0) = area;
    dData(1) = as[0];
    dData(2) = as[1];
    if (theChannel.sendVector(dbTag, commitTag, dData) < 0)  {
	g3ErrorHandler->warning("UniaxialFiber3d::sendSelf() - %s\n",
				"failed to send Vector data");
	return -2;
    }    

    // now invoke sendSelf on the material
    if (theMaterial->sendSelf(commitTag, theChannel) < 0) {
	g3ErrorHandler->warning("UniaxialFiber3d::sendSelf() - %s\n",
				"the material failed in sendSelf()");
	return -3;
    }    	
    
    return 0;
}


int   
UniaxialFiber3d::recvSelf(int commitTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
    // 
    // get tag and material info from an ID
    //

    static ID idData(3);
    int dbTag = this->getDbTag();
    
    if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
	g3ErrorHandler->warning("UniaxialFiber3d::recvSelf() - %s\n",
			      "failed to recv ID data");
	return -1;
    }    

    this->setTag(idData(0));

    // 
    // get area and position datafrom a vector
    //
    
    static Vector dData(3);
    if (theChannel.recvVector(dbTag, commitTag, dData) < 0)  {
	g3ErrorHandler->warning("UniaxialFiber3d::recvSelf() - %s\n",
				"failed to recv Vector data");
	return -2;
    }        
    area = dData(0);
    as[0] = dData(1);
    as[1] = dData(2);

    //
    // now we do the material stuff
    //
    
    int matClassTag = idData(1);    
    
    // if we have a material, check it is of correct type
    if (theMaterial != 0) {
	if (matClassTag != theMaterial->getClassTag()) {
	    delete theMaterial;
	    theMaterial = 0;
	} 
    }

    // if no material we need to get one,
    // NOTE: not an else if in case deleted in if above
    if (theMaterial == 0) {
	theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
	if (theMaterial == 0) {
	    g3ErrorHandler->warning("UniaxialFiber3d::recvSelf() - %s %d\n",
				    "failed to get a UniaxialMaterial of type",
				    matClassTag);
	    return -3;
	}
    }

    // set the materials dbTag and invoke recvSelf on the material
    theMaterial->setDbTag(idData(2));

    // now invoke recvSelf on the material
    if (theMaterial->recvSelf(commitTag, theChannel, theBroker) < 0) {
	g3ErrorHandler->warning("UniaxialFiber3d::recvSelf() - %s\n",
				"the material failed in recvSelf()");
	return -4;
    }    	

    return 0;
}


void UniaxialFiber3d::Print(ostream &s, int flag)
{
    
    if (flag == 1) {
    s << "as[1 2] " << as[0] << " " << as[1] << " Material: " << theMaterial->getTag() << " " << theMaterial->getStrain() << " " << theMaterial->getTangent() << endl;
  }
//    s << "\nFiber: " << this->getTag() << " Type: UniaxialFiber3d";
//    s << "\tMaterial: " << *theMaterial;
//    s << "\tMatrix as:\n" << as; 
//    s << "\tArea: " << area; 

    char st[200];
    sprintf (st, "Material=  %3d  y=%12.4e  z=%12.4e  A=%12.4e\n",theMaterial->getTag(),
             -as[0],as[1],area);   

//    s << " Material: " << theMaterial->getTag();
//    s << "  y: " << -as[1] << " z: " <<  as[2];
//    s << "  Area: " << area << endl; 
      s << st;
}
   
   
ostream &operator<<(ostream &s, UniaxialFiber3d &E)  
{
    E.Print(s); 
    return s; 
}

