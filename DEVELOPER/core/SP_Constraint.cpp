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
                                                                        
// $Revision: 1.6 $
// $Date: 2010-04-23 22:50:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/SP_Constraint.cpp,v $
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the implementation of class SP_Constraint.

#include <SP_Constraint.h>
#include <classTags.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string>
#include <elementAPI.h>
#include <Domain.h>
#include <Node.h>
#include <ID.h>

static int numSPs = 0;
static int nextTag = 0;

int OPS_HomogeneousBC()
{
    Domain* theDomain = OPS_GetDomain();
    if(theDomain == 0) {
	opserr<<"WARNING: domain is not defined\n";
	return -1;
    }
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"insufficient number of args\n";
	return -1;
    }

    // get tag and constr values
    int num = OPS_GetNumRemainingInputArgs();
    ID vals(num);
    if(OPS_GetIntInput(&num, &vals(0)) < 0) {
	opserr << "WARNING invalid int values\n";
	return -1;
    }

    // get node
    Node* theNode = theDomain->getNode(vals(0));
    if(theNode == 0) {
	opserr<<"ERROR node "<<vals(0)<<" is not defined\n";
	return -1;
    }
    int ndf = theNode->getNumberDOF();

    // create homogeneous constraints
    if(vals.Size()-1 < ndf) {
	opserr<<"WARNING: invalid # of constraint values\n";
	return -1;
    }
    for(int i=0; i<ndf; i++) {
	if(vals(i+1) == 0) continue;
	SP_Constraint* theSP = new SP_Constraint(vals(0), i, 0.0, true);
	if(theSP == 0) {
	    opserr<<"WARNING: failed to create SP\n";
	    return -1;
	}
	if(theDomain->addSP_Constraint(theSP) == false) {
	    opserr<<"WARNING: failed to add SP to domain\n";
	    delete theSP;
	    return -1;
	}
    }

    return 0;
}

int OPS_HomogeneousBC_X()
{
    Domain* theDomain = OPS_GetDomain();
    if(theDomain == 0) {
	opserr<<"WARNING: domain is not defined\n";
	return -1;
    }
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"insufficient number of args\n";
	return -1;
    }

    // get the xCrd of nodes to be constrained
    double xLoc;
    int numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &xLoc) < 0) {
	opserr << "WARNING invalid xLoc\n";
	return -1;
    }

    // read in the fixities
    ID fixity(0,3);
    while (OPS_GetNumRemainingInputArgs() > 0) {
	int fix;
	if (OPS_GetIntInput(&numdata, &fix) < 0) {
	    // back one arg
	    OPS_ResetCurrentInputArg(-1);
	    break;
	}
	fixity[fixity.Size()] = fix;
    }

    // set the tolerance, the allowable difference in nodal coordinate and
    // what the value user specified to see if node is constrained or not
    double tol = 1e-10;
    if (OPS_GetNumRemainingInputArgs() > 1) {
	const char* arg = OPS_GetString();
	if (strcmp(arg, "-tol") == 0) {
	    if (OPS_GetDoubleInput(&numdata, &tol) < 0) {
		opserr << "WARNING invalid tol\n";
		return -1;
	    }
	}
    }

    theDomain->addSP_Constraint(0, xLoc, fixity, tol);

    return 0;
}

int OPS_HomogeneousBC_Y()
{
    Domain* theDomain = OPS_GetDomain();
    if(theDomain == 0) {
	opserr<<"WARNING: domain is not defined\n";
	return -1;
    }
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"insufficient number of args\n";
	return -1;
    }

    // get the yCrd of nodes to be constrained
    double yLoc;
    int numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &yLoc) < 0) {
	opserr << "WARNING invalid yLoc\n";
	return -1;
    }

    // read in the fixities
    ID fixity(0,3);
    while (OPS_GetNumRemainingInputArgs() > 0) {
	int fix;
	if (OPS_GetIntInput(&numdata, &fix) < 0) {
	    // back one arg
	    OPS_ResetCurrentInputArg(-1);
	    break;
	}
	fixity[fixity.Size()] = fix;
    }

    // set the tolerance, the allowable difference in nodal coordinate and
    // what the value user specified to see if node is constrained or not
    double tol = 1e-10;
    if (OPS_GetNumRemainingInputArgs() > 1) {
	const char* arg = OPS_GetString();
	if (strcmp(arg, "-tol") == 0) {
	    if (OPS_GetDoubleInput(&numdata, &tol) < 0) {
		opserr << "WARNING invalid tol\n";
		return -1;
	    }
	}
    }

    theDomain->addSP_Constraint(1, yLoc, fixity, tol);

    return 0;
}

int OPS_HomogeneousBC_Z()
{
    Domain* theDomain = OPS_GetDomain();
    if(theDomain == 0) {
	opserr<<"WARNING: domain is not defined\n";
	return -1;
    }
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"insufficient number of args\n";
	return -1;
    }

    // get the zCrd of nodes to be constrained
    double zLoc;
    int numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &zLoc) < 0) {
	opserr << "WARNING invalid zLoc\n";
	return -1;
    }

    // read in the fixities
    ID fixity(0,3);
    while (OPS_GetNumRemainingInputArgs() > 0) {
	int fix;
	if (OPS_GetIntInput(&numdata, &fix) < 0) {
	    // back one arg
	    OPS_ResetCurrentInputArg(-1);
	    break;
	}
	fixity[fixity.Size()] = fix;
    }

    // set the tolerance, the allowable difference in nodal coordinate and
    // what the value user specified to see if node is constrained or not
    double tol = 1e-10;
    if (OPS_GetNumRemainingInputArgs() > 1) {
	const char* arg = OPS_GetString();
	if (strcmp(arg, "-tol") == 0) {
	    if (OPS_GetDoubleInput(&numdata, &tol) < 0) {
		opserr << "WARNING invalid tol\n";
		return -1;
	    }
	}
    }

    theDomain->addSP_Constraint(2, zLoc, fixity, tol);

    return 0;
}

// 2 little procedures needed for parallel processing all due to fact that SP's need 
// to keep unique tags among processes in parallel

int SP_Constraint_SetNextTag(int next) {
  nextTag = next;
  return nextTag;
}

int SP_Constraint_GetNextTag(void) {
  return nextTag;
}

// constructor for FEM_ObjectBroker
SP_Constraint::SP_Constraint(int clasTag)
:DomainComponent(0,clasTag),
 nodeTag(0), dofNumber(0), valueR(0.0), valueC(0.0), isConstant(true), 
 loadPatternTag(-1)
{
  numSPs++;
}

// constructor for a subclass to use
SP_Constraint::SP_Constraint(int node, int ndof, int clasTag)
:DomainComponent(nextTag++, clasTag),
 nodeTag(node), dofNumber(ndof), valueR(0.0), valueC(0.0), isConstant(true), 
 loadPatternTag(-1)
 // valueC is set to 1.0 so that homo will be false when recvSelf() invoked
 // should be ok as valueC cannot be used by subclasses and subclasses should
 // not be used if it is a homogeneous constraint.
{
  numSPs++;
}

// constructor for object of type SP_Constraint
SP_Constraint::SP_Constraint(int node, int ndof, double value, bool ISconstant)
:DomainComponent(nextTag++, CNSTRNT_TAG_SP_Constraint),
 nodeTag(node), dofNumber(ndof), valueR(value), valueC(value), isConstant(ISconstant),
 loadPatternTag(-1)
{
  numSPs++;
}

SP_Constraint::~SP_Constraint()
{
  numSPs--;
  if (numSPs == 0)
    nextTag = 0;
}

int
SP_Constraint::getNodeTag(void) const
{
    // return id of constrained node
    return nodeTag;
}

int
SP_Constraint::getDOF_Number(void) const
{
    //  return the number of the constrained DOF    
    return dofNumber;
}


double
SP_Constraint::getValue(void)
{
    // return the value of the constraint
    return valueC;
}

int
SP_Constraint::applyConstraint(double loadFactor)
{
    // as SP_Constraint objects are time invariant nothing is done
    if (isConstant == false)
	valueC = loadFactor*valueR;

    return 0;
}


bool
SP_Constraint::isHomogeneous(void) const
{
    if (valueR == 0.0)
	return true;
    else
	return false;
}

void
SP_Constraint::setLoadPatternTag(int tag)
{
  loadPatternTag = tag;
}

int
SP_Constraint::getLoadPatternTag(void) const
{
  return loadPatternTag;
}

int 
SP_Constraint::sendSelf(int cTag, Channel &theChannel)
{
    static Vector data(8);  // we send as double to avoid having 
                     // to send two messages.
    data(0) = this->getTag(); 
    data(1) = nodeTag;
    data(2) = dofNumber;
    data(3) = valueC;
    if (isConstant == true)
	data(4) = 1.0;
    else
	data(4) = 0.0;
    data(5) = valueR;
    data(6) = this->getLoadPatternTag();

    data(7) = nextTag;

    int result = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (result != 0) {
      opserr << "WARNING SP_Constraint::sendSelf - error sending Vector data\n";
      return result;
    }

    return 0;
}

int 
SP_Constraint::recvSelf(int cTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
    static Vector data(8);  // we sent the data as double to avoid having to send
                     // two messages
    int result = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (result < 0) {
	opserr << "WARNING SP_Constraint::recvSelf - error receiving Vector data\n";
	return result;
    }
    
    // if o.k. set the data
    this->setTag((int)data(0));
    nodeTag = (int)data(1);
    dofNumber = (int)data(2);
    valueC = data(3);

    if (data(4) == 1.0)
	isConstant = true;
    else
	isConstant = false;
    valueR = data(5);
    valueC = valueR;
    this->setLoadPatternTag((int)data(6));

    nextTag = (int)data(7);

    return 0;
}


void
SP_Constraint::Print(OPS_Stream &s, int flag) 
{
    s << "SP_Constraint: " << this->getTag();
    s << "\t Node: " << nodeTag << " DOF: " << dofNumber+1;
    s << " ref value: " << valueR << " current value: " << valueC << endln;
}








