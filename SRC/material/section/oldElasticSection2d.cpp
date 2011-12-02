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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/oldElasticSection2d.cpp,v $
                                                                        
                                                                        
///////////////////////////////////////////////////////
// File:  ~/Src/element/hinge/ElasticSection2d.cpp
//
// Written by Matthew Peavy
//
// Written:	 Feb 13, 2000
// Debugged: Feb 14, 2000
// Revised:		   , 200x
//
//
// Purpose:  This file contains the function definitions
// for the ElasticSection2d class.


#include <ElasticSection2d.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <stdlib.h>

ID ElasticSection2d::code(2);

ElasticSection2d::ElasticSection2d(void)
:SectionForceDeformation(0, SECTION_TAG_ElasticSection2d),
 E(0), I(0), A(0),
 k(2,2), f(2,2),
 e(2), s(2),
 massDens(0)
{
}

ElasticSection2d::ElasticSection2d (int tag, 
				    double E_in, double I_in, double A_in,
				    double massDensityPerUnitLength)	  //massDens defaults to 0 in .h file
:SectionForceDeformation(tag, SECTION_TAG_ElasticSection2d),
 E(E_in), I(I_in), A(A_in),
 k(2,2), f(2,2),
 e(2), s(2),
 massDens(massDensityPerUnitLength)
{	
    if (E <= 0.0)  {
	cerr << "FATAL - ElasticSection2d::ElasticSection2d() - Paramater E is zero or negative.";
	//exit(-1); 
    }
    
    if (I <= 0.0)  {
	cerr << "FATAL - ElasticSection2d::ElasticSection2d() - Parameter I is zero or negative.";
	//exit(-1); 
    }
    
    if (A <= 0.0)  {
	cerr << "FATAL - ElasticSection2d::ElasticSection2d() - Parameter A is zero or negative.";
	//exit(-1); 
    }

    k(0,0) = E*A;
    k(1,1) = E*I;
    
    f(0,0) = 1/(E*A);
    f(1,1) = 1/(E*I);
    
    if (code(0) == 0)
    {
	code(0) = 2;	// Px is the first quantity
	code(1) = 1;	// Mz is the second
    }
}	

ElasticSection2d::~ElasticSection2d(void)
{
} 

int 
ElasticSection2d::commitState(void)
{
    return 0;	// Elastic Hinge -- nothing to commit
}

int 
ElasticSection2d::revertToLastCommit(void)
{
    return 0;	//Elastic Hinge -- nothing to revert to
}

void
ElasticSection2d::setTrialDeformation (const Vector &def)
{
    e = def;
    return;
}

const Vector &
ElasticSection2d::getDeformation (void)
{
    return e;
}

const Vector &
ElasticSection2d::getResistingForce (void)
{
    s = k*e;
    return s;
}

const Vector &
ElasticSection2d::getPrevResistingForce (void)
{
    return s;
}

const Matrix &
ElasticSection2d::getTangentStiff(void)
{
    return k;
}

const Matrix &
ElasticSection2d::getPrevTangentStiff (void)
{
    return k;
}

const Matrix &
ElasticSection2d::getFlexMatrix (void)
{
    return f;
}

const Matrix &
ElasticSection2d::getPrevFlexMatrix (void)
{
    return f;
}

SectionForceDeformation*
ElasticSection2d::getCopy ()
{
    // Make a copy of the hinge
    ElasticSection2d *theCopy = new ElasticSection2d (this->getTag(), E, I, A, massDens);

    theCopy->e = e;
    theCopy->s = s;
    
    return theCopy;
}

const ID&
ElasticSection2d::getType () const
{
    return code;
}

int
ElasticSection2d::getOrder () const
{
    return 2;
}

int
ElasticSection2d::sendSelf(int commitTag, Channel &theChannel)
{
    Vector data(5);
    int result = 0;
    int dataTag = this->getDbTag();
    
    data(0) = E;
    data(1) = I;
    data(2) = A;
    data(3) = massDens;
    data(4) = this->getTag();
    
    result = theChannel.sendVector(dataTag, commitTag, data);
    if(result<0) {
	cerr << "ElasticSection2d::sendSelf - Failed to send data.\n";
	return -1;
    }
    
    return 0;
}

int
ElasticSection2d::recvSelf(int commitTag, Channel &theChannel,
					 FEM_ObjectBroker &theBroker)
{
    Vector data(5);
    int result = 0;
    int dataTag = this->getDbTag();
    
    result = theChannel.recvVector(dataTag, commitTag, data);
    if(result <0) {
	cerr << "ElasticSection2d::recvSelf - Failed to receive data.\n";
	return -1;
    }
    
    E = data(0);
    I = data(1);
    A = data(2);
    massDens = data(3);
    this->setTag((int)data(4));
    
    return 0;
}
 
void
ElasticSection2d::Print(ostream &s, int flag)
{
    s << "\nHinge tag: " <<this->getTag() << "   Type: ElasticSection2d\n";
    s << "Properties:   E = " << E << endl;
    s << "\tI = " << I << endl;
    s << "\tA = " << A << endl;
}

ostream &operator<<(ostream &s, ElasticSection2d &H)
{
	H.Print(s);
	return s;
}
