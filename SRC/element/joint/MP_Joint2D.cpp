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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-04-02 22:02:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/MP_Joint2D.cpp,v $

// Written: Arash & GGD
// Created: 08/01
// Revision: Arash

// Purpose: This file contains the implementation of class MP_TimeVary.


#include <MP_Joint2D.h>

#include <stdlib.h>
#include <math.h>

#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
 
				// main degree of freedom for rotation

// constructor for FEM_ObjectBroker
MP_Joint2D::MP_Joint2D()
:MP_Constraint( 0 , CNSTRNT_TAG_MP_Joint2D ),thisDomain(0),
 nodeRetained(0),nodeConstrained(0), MainDOF(0), AuxDOF(0), FixedEnd(0),
 constraint(0), constrDOF(0),retainDOF(0),dbTag1(0), dbTag2(0), dbTag3(0),
 RetainedNode(0), ConstrainedNode(0), LargeDisplacement(0), Length0(0.0)
{
    
}


// general constructor for ModelBuilder
MP_Joint2D::MP_Joint2D(Domain *theDomain, int tag, int nodeRetain, int nodeConstr,
		int Maindof, int fixedend, int LrgDsp )
:MP_Constraint( tag , CNSTRNT_TAG_MP_Joint2D ), thisDomain(theDomain),
 nodeRetained(nodeRetain), nodeConstrained(nodeConstr), MainDOF(Maindof), AuxDOF(0),
 FixedEnd(fixedend), constraint(0), constrDOF(0), retainDOF(0),
 dbTag1(0), dbTag2(0), dbTag3(0), RetainedNode(0), ConstrainedNode(0),
 LargeDisplacement( LrgDsp ), Length0(0.0)
{

  if( thisDomain == NULL ) {
    opserr << "WARNING MP_Joint2D(): Specified domain does not exist";
    opserr << "Domain = 0\n";
    return;
  }

  this->setTag(tag);

	// get node pointers of constrainted and retained nodes
	ConstrainedNode = theDomain->getNode(nodeConstrained);
	if (ConstrainedNode == NULL)
	{
		opserr << "MP_Joint2D::MP_Joint2D: nodeConstrained: ";
		opserr << nodeConstrained << "does not exist in model\n";
		exit(0);
	}

	RetainedNode = theDomain->getNode(nodeRetained);
	if (RetainedNode == NULL)
	{
		opserr << "MP_Joint2D::MP_Joint2D: nodeRetained: ";
		opserr << nodeRetained << "does not exist in model\n";
		exit(0);
	}

	// check for proper degrees of freedom
	int RnumDOF = RetainedNode->getNumberDOF();
	int CnumDOF = ConstrainedNode->getNumberDOF();
    if (RnumDOF != 4 || CnumDOF != 3 ){
		opserr << "MP_Joint2D::MP_Joint2D - mismatch in numDOF\n DOF not supported by this type of constraint";
		return;
    }

	// check the main degree of freedom. Assign auxilary DOF 
	if ( MainDOF!= 2 && MainDOF!=3 ) {
			opserr << "MP_Joint2D::MP_Joint2D - Wrong main degree of freedom" ;
			return;
    }
	if ( MainDOF == 2 ) AuxDOF = 3;
	if ( MainDOF == 3 ) AuxDOF = 2;
	
	// check the fixed end flag
	if ( FixedEnd!= 0 && FixedEnd!=1 ) {
			opserr << "MP_Joint2D::MP_Joint2D - Wrong fixed end flag";
			return;
    }
	

	// check for proper dimensions of coordinate space
	const Vector &crdR = RetainedNode->getCrds();
    int dimR = crdR.Size();
	const Vector &crdC = ConstrainedNode->getCrds();
    int dimC = crdC.Size();
    
	if (dimR != 2 || dimC != 2 ){
		opserr << "MP_Joint2D::MP_Joint2D - mismatch in dimnesion\n dimension not supported by this type of constraint";
		return;
    }


	// calculate the initial length of the rigid link
	double deltaX = crdC(0) - crdR(0);
	double deltaY = crdC(1) - crdR(1);

	Length0 = sqrt( deltaX*deltaX + deltaY*deltaY );
    if ( Length0 <= 1.0e-12 ) {
		opserr << "MP_Joint2D::MP_Joint2D - The constraint length is zero\n";
    }
   
	// allocate and set up the constranted and retained id's
	// allocate and set up the constraint matrix
	if ( FixedEnd == 0 )
	{
		// the end is released
		constrDOF = new ID(CnumDOF-1);
		retainDOF = new ID(RnumDOF-1);
		
		(*constrDOF)(0) = 0;
		(*constrDOF)(1) = 1;

		(*retainDOF)(0) = 0;
		(*retainDOF)(1) = 1;
		(*retainDOF)(2) = MainDOF;
		
		constraint = new Matrix( CnumDOF-1 , RnumDOF-1 );
		
		(*constraint) (0,0) = 1.0 ;
		(*constraint) (0,2) = -deltaY ;
		(*constraint) (1,1) = 1.0 ;
		(*constraint) (1,2) = deltaX ;
	} else
	{
		// the end is fixed
		constrDOF = new ID(CnumDOF);
		retainDOF = new ID(RnumDOF);
		
		(*constrDOF)(0) = 0;
		(*constrDOF)(1) = 1;
		(*constrDOF)(2) = 2;
		
		(*retainDOF)(0) = 0;
		(*retainDOF)(1) = 1;
		(*retainDOF)(2) = 2;
		(*retainDOF)(3) = 3;
		
		constraint = new Matrix( CnumDOF , RnumDOF );
		
		(*constraint) (0,0) = 1.0 ;
		(*constraint) (0,MainDOF) = -deltaY ;
		(*constraint) (1,1) = 1.0 ;
		(*constraint) (1,MainDOF) = deltaX ;
		(*constraint) (2,AuxDOF) = 1.0 ;
	}
 
	if (constrDOF == NULL || retainDOF == NULL ) { 
		opserr << "MP_Joint2D::MP_Joint2D - ran out of memory \ncan not generate ID for nodes\n";
		exit(-1);
	}
	
	if (constraint == NULL ) {
		opserr << "MP_Joint2D::MP_Joint2D - ran out of memory \ncan not generate the constraint matrix";
		exit(-1);
    }
}



MP_Joint2D::~MP_Joint2D()
{
    // invoke the destructor on the matrix and the two ID objects
    if (constraint != NULL)
	delete constraint;
    if (constrDOF != NULL)
	delete constrDOF;
    if (retainDOF != NULL)
	delete retainDOF;    
}


int
MP_Joint2D::getNodeRetained(void) const
{
    // return id of retained node
    return nodeRetained;
}

int
MP_Joint2D::getNodeConstrained(void) const
{
    // return id of constrained node    
    return nodeConstrained;
}


const ID &
MP_Joint2D::getConstrainedDOFs(void) const
{
    if (constrDOF == NULL) {
	opserr << "MP_Joint2D::getConstrainedDOF - no ID was set, ";
	opserr << "was recvSelf() ever called? or subclass incorrect?\n";	
	exit(-1);
    }

    // return the ID corresponding to constrained DOF of Ccr
    return (*constrDOF);    
}


const ID &
MP_Joint2D::getRetainedDOFs(void) const
{
    if (retainDOF == NULL) {
	opserr << "MP_Joint2D::getRetainedDOFs - no ID was set\n ";
	opserr << "was recvSelf() ever called? or subclass incorrect?\n";		
	exit(-1);
    }

    // return the ID corresponding to retained DOF of Ccr
    return (*retainDOF);    
}


int 
MP_Joint2D::applyConstraint(double timeStamp)
{
    if ( LargeDisplacement != 0 )
	{
		// calculate the constraint at this moment

		// get the coordinates of the two nodes - check dimensions are the same FOR THE MOMENT
		const Vector &crdR = RetainedNode->getCrds();
		const Vector &crdC = ConstrainedNode->getCrds();

		// get commited displacements of nodes to get updated coordinates
		const Vector &dispR = RetainedNode->getDisp();
		const Vector &dispC = ConstrainedNode->getDisp();

		double deltaX = dispC(0) + crdC(0) - dispR(0) - crdR(0);
		double deltaY = dispC(1) + crdC(1) - dispR(1) - crdR(1);

		constraint->Zero();
		
		if ( FixedEnd == 0 )
		{
			// the end is released
			(*constraint) (0,0) = 1.0 ;
			(*constraint) (0,2) = -deltaY ;
			(*constraint) (1,1) = 1.0 ;
			(*constraint) (1,2) = deltaX ;
		} else
		{
			// the end is fixed
			(*constraint) (0,0) = 1.0 ;
			(*constraint) (0,MainDOF) = -deltaY ;
			(*constraint) (1,1) = 1.0 ;
			(*constraint) (1,MainDOF) = deltaX ;
			(*constraint) (2,AuxDOF) = 1.0 ;
		}
	}
	return 0;
}


bool
MP_Joint2D::isTimeVarying(void) const
{
    if ( LargeDisplacement != 0 ) return true;

	return false;
}


int MP_Joint2D::sendSelf(int commitTag, Channel &theChannel)
{
	Vector data(15);
    int dataTag = this->getDbTag();

    data(0) = this->getTag(); 
    data(1) = nodeRetained;
    data(2) = nodeConstrained;
    data(3) = MainDOF;
	data(4) = AuxDOF;
	data(5) = FixedEnd;
    
	if (constrDOF == 0) data(6) = 0; else data(6) = constrDOF->Size();    
	if (retainDOF == 0) data(7) = 0; else data(7) = retainDOF->Size();        
    if (constraint == 0) data(8) = 0; else data(8) = constraint->noRows();
	if (constraint == 0) data(9) = 0; else data(9) = constraint->noCols();   
    // need two database tags for ID objects
    if (constrDOF != 0 && dbTag1 == 0) dbTag1 = theChannel.getDbTag();
    if (retainDOF != 0 && dbTag2 == 0) dbTag2 = theChannel.getDbTag();
	if (constraint != 0 && dbTag3 == 0) dbTag3 = theChannel.getDbTag();

    data(10) = dbTag1;
    data(11) = dbTag2;
	data(12) = dbTag3;
    data(13) = LargeDisplacement;
    data(14) = Length0;

	// now send the data vector
    int result = theChannel.sendVector(dataTag, commitTag, data);
    if (result < 0) {
		opserr << "WARNING MP_Joint2D::sendSelf - error sending ID data\n";
		return result;  
    }    
    
	// send constrDOF
    if (constrDOF != 0 && constrDOF->Size() != 0) {
		int result = theChannel.sendID(dbTag1, commitTag, *constrDOF);
		if (result < 0) {
			opserr << "WARNING MP_Joint2D::sendSelf ";
			opserr << "- error sending constrained DOF data\n";
			return result;
		}
	}

	// send retainDOF
    if (retainDOF != 0 && retainDOF->Size() != 0) {
		int result = theChannel.sendID(dbTag2, commitTag, *retainDOF);
		if (result < 0) {
			opserr << "WARNING MP_Joint2D::sendSelf ";
			opserr << "- error sending retained DOF data\n";
			return result;
		}
    }

	// send constraint matrix 
    if (constraint != 0 && constraint->noRows() != 0) {


	int result = theChannel.sendMatrix(dbTag3, commitTag, *constraint);
	if (result < 0) {
	    opserr << "WARNING MP_Joint2D::sendSelf ";
	    opserr << "- error sending constraint Matrix data\n"; 
	    return result;  
	}
    }

    return 0;
}


int MP_Joint2D::recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    Vector data(15);
    int result = theChannel.recvVector(dataTag, commitTag, data);
    if (result < 0) {
	opserr << "WARNING MP_Joint2D::recvSelf - error receiving ID data\n";
	return result;  
    }    

    this->setTag( (int) data(0));
    nodeRetained = (int) data(1);
    nodeConstrained = (int) data(2);
    MainDOF = (int) data(3);
	AuxDOF = (int) data(4);
	FixedEnd = (int) data(5);

	int constrDOFsize = (int) data(6);
	int retainDOFsize = (int) data(7);
    int numRows = (int) data(8);
	int numCols = (int) data(9); 
	
	dbTag1 = (int) data(10);
    dbTag2 = (int) data(11);
    dbTag3 = (int) data(12);
	LargeDisplacement = (int) data(13);
    Length0 = data(14);


    // receive constrDOF ID
    if (constrDOFsize != 0) {
	constrDOF = new ID(constrDOFsize);
	int result = theChannel.recvID(dbTag1, commitTag, *constrDOF);
	if (result < 0) {
	    opserr << "WARNING MP_Joint2D::recvSelf ";
	    opserr << "- error receiving constrained data\n"; 
	    return result;  
	}	
    }
    
    // receive retainDOF ID
    if (retainDOFsize != 0) {
	retainDOF = new ID(retainDOFsize);
	int result = theChannel.recvID(dbTag2, commitTag, *retainDOF);
	if (result < 0) {
	    opserr << "WARNING MP_Joint2D::recvSelf ";
	    opserr << "- error receiving retained data\n"; 
	    return result;  
	}	
    }    
    
    // receive the constraint matrix
	if (numRows != 0 && numCols != 0) {
		constraint = new Matrix(numRows,numCols);

		int result = theChannel.recvMatrix(dbTag3, commitTag, *constraint);
		
		if (result < 0) {
			opserr << "WARNING MP_Joint2D::recvSelf ";
			opserr << "- error receiving Matrix data\n";
			return result;
		}
	}

    return 0;
}


const Matrix &MP_Joint2D::getConstraint(void)
{
    if (constraint == 0) {
	opserr << "MP_Joint2D::getConstraint - no Matrix was set\n";
	exit(-1);
    }    

	// Length correction
	// to correct the trial displacement
    if ( LargeDisplacement == 2 )
	{
		// get the coordinates of the two nodes - check dimensions are the same FOR THE MOMENT
		const Vector &crdR = RetainedNode->getCrds();
		const Vector &crdC = ConstrainedNode->getCrds();

		// get commited displacements of nodes to get updated coordinates
		const Vector &dispR = RetainedNode->getTrialDisp();
		const Vector &dispC = ConstrainedNode->getTrialDisp();

		double deltaX = dispC(0) + crdC(0) - dispR(0) - crdR(0);
		double deltaY = dispC(1) + crdC(1) - dispR(1) - crdR(1);


		Vector Direction(2);
		Direction(0) = deltaX;
		Direction(1) = deltaY;
		double NewLength = Direction.Norm();
		if ( NewLength < 1e-12 ) opserr << "MP_Joint2D::applyConstraint : length of rigid link is too small or zero"; 
		Direction = Direction * (Length0/NewLength);		// correct the length
		// find new displacements of the constrainted node
	
		Vector NewLocation(3);
		NewLocation(0) = Direction(0) + dispR(0) + crdR(0) - crdC(0);
		NewLocation(1) = Direction(1) + dispR(1) + crdR(1) - crdC(1);
		NewLocation(2) = dispC(2);
		int dummy = ConstrainedNode->setTrialDisp( NewLocation );
	}
	// end of length correction procedure

    // return the constraint matrix Ccr
    return (*constraint);
}
    
void MP_Joint2D::Print(OPS_Stream &s, int flag )
{
    s << "MP_Joint2D: " << this->getTag() << "\n";
    s << "\tConstrained Node: " << nodeConstrained;
    s << " Retained Node: " << nodeRetained ;
	s << " Fixed end: " << FixedEnd << " Large Disp: " << LargeDisplacement;
    if (constrDOF != 0)
	s << " constrained dof: " << *constrDOF;    
    if (retainDOF != 0)
	s << " retained dof: " << *retainDOF;        
    if (constraint != 0)
	s << " constraint matrix: " << *constraint << "\n";

}


void
MP_Joint2D::setDomain(Domain *theDomain)
{
	this->DomainComponent::setDomain(theDomain);
	thisDomain = theDomain;

    RetainedNode = thisDomain->getNode(nodeRetained);
    ConstrainedNode = thisDomain->getNode(nodeConstrained);
}


