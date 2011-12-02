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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/oldElasticSection3d.cpp,v $
                                                                        
                                                                        
///////////////////////////////////////////////////////
// File:  ~/Src/element/hinge/ElasticSection3d.cpp
//
// Written by Matthew Peavy
//
// Written:	 Feb 13, 2000
// Debugged: Feb 14, 2000
// Revised:		   , 200x
//
//
// Purpose:  This file contains the function definitions
// for the ElasticSection3d class.  

#include <ElasticSection3d.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <stdlib.h>

ID ElasticSection3d::code(6);

ElasticSection3d::ElasticSection3d(void)
:SectionForceDeformation(0, SECTION_TAG_ElasticSection3d),
 E(0), Iz(0), Iy(0), A(0),
 k(3,3), f(3,3),
 e(3), s(3),
 massDens(0)
{
}

ElasticSection3d::ElasticSection3d (int tag, 
				    double E_in, double Iz_in, double Iy_in,
				    double A_in,
				    double massDensityPerUnitLength)	  //massDens defaults to 0 in .h file
:SectionForceDeformation(tag, SECTION_TAG_ElasticSection3d),
 E(E_in), Iz(Iz_in), Iy(Iy_in), A(A_in),
 k(3,3), f(3,3),
 e(3), s(3),
 massDens(massDensityPerUnitLength)
{	
    if (E <= 0.0)  {
	cerr << "FATAL - ElasticSection3d::ElasticSection3d() - Paramater E is zero or negative.";
	//exit(-1); 
    }
    
    if (Iz <= 0.0)  {
	cerr << "FATAL - ElasticSection3d::ElasticSection3d() - Parameter Iz is zero or negative.";
	//exit(-1); 
    }
    
    if (Iy <= 0.0)  {
	cerr << "FATAL - ElasticSection3d::ElasticSection3d() - Parameter Iy is zero or negative.";
	//exit(-1); 
    }
    
    if (A <= 0.0)  {
	cerr << "FATAL - ElasticSection3d::ElasticSection3d() - Parameter A is zero or negative.";
	//exit(-1); 
    }

    k(0,0) = E*A;
    k(1,1) = E*Iz;
    k(2,2) = E*Iy;
    
    f(0,0) = 1/(E*A);
    f(1,1) = 1/(E*Iz);
    f(2,2) = 1/(E*Iy);
    
    if (code(0) == 0)
    {
	code(0) = 2;	// Px is the first quantity
	code(1) = 1;	// Mz is the second
	code(2) = 4;	// My is the third
    }
}	

ElasticSection3d::~ElasticSection3d(void)
{
} 

int 
ElasticSection3d::commitState(void)
{
    return 0;	// Elastic Hinge -- nothing to commit
}

int 
ElasticSection3d::revertToLastCommit(void)
{
    return 0;	//Elastic Hinge -- nothing to revert to
}

void
ElasticSection3d::setTrialDeformation (const Vector &def)
{
    e = def;
    return;
}

const Vector &
ElasticSection3d::getDeformation (void)
{
    return e;
}

const Vector &
ElasticSection3d::getResistingForce (void)
{
    s = k*e;
    return s;
}

const Vector &
ElasticSection3d::getPrevResistingForce (void)
{
    return s;
}

const Matrix &
ElasticSection3d::getTangentStiff(void)
{
    return k;
}

const Matrix &
ElasticSection3d::getPrevTangentStiff (void)
{
    return k;
}

const Matrix &
ElasticSection3d::getFlexMatrix (void)
{
    return f;
}

const Matrix &
ElasticSection3d::getPrevFlexMatrix (void)
{
    return f;
}

SectionForceDeformation*
ElasticSection3d::getCopy ()
{
    // Make a copy of the hinge
    ElasticSection3d *theCopy =
	new ElasticSection3d (this->getTag(), E, Iz, Iy, A, massDens);

    theCopy->e = e;
    theCopy->s = s;
    
    return theCopy;
}

const ID&
ElasticSection3d::getType () const
{
    return code;
}

int
ElasticSection3d::getOrder () const
{
    return 3;
}

int
ElasticSection3d::sendSelf(int commitTag, Channel &theChannel)
{
    Vector data(6);
    int result = 0;
    int dataTag = this->getDbTag();
    
    data(0) = E;
    data(1) = Iz;
    data(2) = Iy;
    data(3) = A;
    data(4) = massDens;
    data(5) = this->getTag();
    
    result = theChannel.sendVector(dataTag, commitTag, data);
    if(result<0) {
	cerr << "ElasticSection3d::sendSelf - Failed to send data.\n";
	return -1;
    }
    
    return 0;
}

int
ElasticSection3d::recvSelf(int commitTag, Channel &theChannel,
					 FEM_ObjectBroker &theBroker)
{
    Vector data(6);
    int result = 0;
    int dataTag = this->getDbTag();
    
    result = theChannel.recvVector(dataTag, commitTag, data);
    if(result <0) {
	cerr << "ElasticSection3d::recvSelf - Failed to receive data.\n";
	return -1;
    }
    
    E = data(0);
    Iz = data(1);
    Iy = data(2);    
    A = data(3);
    massDens = data(4);
    this->setTag((int)data(5));
    
    return 0;
}
 
void
ElasticSection3d::Print(ostream &s, int flag)
{
    s << "\nHinge tag: " <<this->getTag() << "   Type: ElasticSection3d\n";
    s << "Properties:   E = " << E << endl;
    s << "\tIz = " << Iz << endl;
    s << "\tIy = " << Iy << endl;
    s << "\tA = " << A << endl;
}

ostream &operator<<(ostream &s, ElasticSection3d &H)
{
    H.Print(s);
    return s;
}
