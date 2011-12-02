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
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticSection2d.cpp,v $
                                                                        
                                                                        
///////////////////////////////////////////////////////
// File:  ~/Src/element/hinge/ElasticSection2d.cpp
//
// Written by Matthew Peavy
//
// Written:  Feb 13, 2000
// Debugged: Feb 14, 2000
// Revised:  May 2000 -- MHS (to include elastic shear)
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

#include <G3Globals.h>
#include <classTags.h>

Vector ElasticSection2d::s(2);
ID ElasticSection2d::code(2);

ElasticSection2d::ElasticSection2d(void)
:SectionForceDeformation(0, SEC_TAG_Elastic2d),
 E(0), A(0), I(0),
 k(2,2), f(2,2), e(2), eCommit(2)
{
    if (code(0) != SECTION_RESPONSE_P)
    {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
    }    
}

ElasticSection2d::ElasticSection2d
(int tag, double E_in, double A_in, double I_in)
:SectionForceDeformation(tag, SEC_TAG_Elastic2d),
 E(E_in), A(A_in), I(I_in),
 k(2,2), f(2,2), e(2), eCommit(2)
{
    if (E <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input E <= 0.0 ... setting E to 1.0",
			"ElasticSection2d::ElasticSection2d");
		E = 1.0;
    }
	
    if (A <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input A <= 0.0 ... setting A to 1.0",
			"ElasticSection2d::ElasticSection2d");
		A = 1.0;
    }
    
    if (I <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input I <= 0.0 ... setting I to 1.0",
			"ElasticSection2d::ElasticSection2d");
		I = 1.0;
    }    
	
    k(0,0) = E*A;
    k(1,1) = E*I;
    
    f(0,0) = 1/k(0,0);
    f(1,1) = 1/k(1,1);

    if (code(0) != SECTION_RESPONSE_P)
    {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
    }
}

ElasticSection2d::ElasticSection2d
(int tag, double EA_in, double EI_in)
:SectionForceDeformation(tag, SEC_TAG_Elastic2d),
 k(2,2), f(2,2), e(2), eCommit(2),
 E(1), A(EA_in), I(EI_in)
{
    if (A <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input EA <= 0.0 ... setting EA to 1.0",
			"ElasticSection2d::ElasticSection2d");
		A = 1.0;
    }
	
    if (I <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input EI <= 0.0 ... setting EI to 1.0",
			"ElasticSection2d::ElasticSection2d");
		I = 1.0;
    }
	
    k(0,0) = A;
    k(1,1) = I;
    
    f(0,0) = 1/A;
    f(1,1) = 1/I;

    if (code(0) != SECTION_RESPONSE_P)
    {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
    }
}

ElasticSection2d::~ElasticSection2d(void)
{

}

int 
ElasticSection2d::commitState(void)
{
	eCommit = e;

    return 0;
}

int 
ElasticSection2d::revertToLastCommit(void)
{
	e = eCommit;

    return 0;
}

int 
ElasticSection2d::revertToStart(void)
{
	eCommit.Zero();

    return 0;
}

int
ElasticSection2d::setTrialSectionDeformation (const Vector &def)
{
    e = def;

    return 0;
}

const Vector &
ElasticSection2d::getSectionDeformation (void)
{
    return e;
}

const Vector &
ElasticSection2d::getStressResultant (void)
{
    // s = k*e;
    s(0) = k(0,0)*e(0);
    s(1) = k(1,1)*e(1);    

    return s;
}

const Matrix &
ElasticSection2d::getSectionTangent(void)
{
    return k;
}

const Matrix &
ElasticSection2d::getSectionFlexibility (void)
{
    return f;
}

SectionForceDeformation*
ElasticSection2d::getCopy ()
{
    // Make a copy of the hinge
    ElasticSection2d *theCopy =
	new ElasticSection2d (this->getTag(), E, A, I);

    theCopy->eCommit = eCommit;

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
	int res = 0;

    static Vector data(6);
    
    int dataTag = this->getDbTag();
    
	data(0) = this->getTag();
    data(1) = E;
    data(2) = A;
    data(3) = I;    
    data(4) = eCommit(0);
	data(5) = eCommit(1);

    res += theChannel.sendVector(dataTag, commitTag, data);
    if (res<0) {
		g3ErrorHandler->warning("%s -- failed to send data",
			"ElasticSection2d::sendSelf");
		return res;
    }
    
    return res;
}

int
ElasticSection2d::recvSelf(int commitTag, Channel &theChannel,
					 FEM_ObjectBroker &theBroker)
{
	int res = 0;

    static Vector data(6);

    int dataTag = this->getDbTag();

    res += theChannel.recvVector(dataTag, commitTag, data);
    if(res < 0) {
		g3ErrorHandler->warning("%s -- failed to receive data",
			"ElasticSection2d::recvSelf");
		return res;
    }
    
	this->setTag((int)data(0));
    E = data(1);
    A = data(2);
    I = data(3);
    eCommit(0) = data(4);
	eCommit(1) = data(5);

	// Set section stiffness and flexibility matrices
	k(0,0) = E*A;
	k(1,1) = E*I;

	f(0,0) = 1/k(0,0);
	f(1,1) = 1/k(1,1);

    return res;
}
 
void
ElasticSection2d::Print(ostream &s, int flag)
{
    s << "ElasticSection2d, tag: " << this->getTag() << endl;
    s << "\tE: " << E << endl;
	s << "\tA: " << A << endl;
	s << "\tI: " << I << endl;
}

