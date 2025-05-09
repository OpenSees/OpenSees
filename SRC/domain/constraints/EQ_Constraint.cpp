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
                                                                        
// $Revision: 1.0 $
// $Date: 2025-05-29$
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/EQ_Constraint.cpp,v $
                                                                        
                                                                        
// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 05/2020
// Revision: A
//
// Purpose: This file contains the implementation of class EQ_Constraint.
//
// The class EQ_Constraint interface:
//

#include <EQ_Constraint.h>

#include <stdlib.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Domain.h>
#include <Node.h>

static int numEQs = 0;
static int nextTag = 0;

int OPS_EquationConstraint()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) {
	    opserr<<"WARNING: domain is not defined\n";
	    return -1;
    }

    int numRemainingArgs = OPS_GetNumRemainingInputArgs();
    if(numRemainingArgs < 6 || numRemainingArgs % 3) {
	    opserr<<"WARNING: invalid # of args: equationConstraint cNodeTag cDOF cCoef rNodeTag1 rDOF1 rCoef1 rNodeTag2 rDOF2 rCoef2 ...\n";
	    return -1;
    }

    int numData = 1;
    int cNode, cDOF;
    if(OPS_GetIntInput(&numData, &cNode)) {
        opserr<<"WARNING invalid cNodeTag inputs\n";
        return -1;
    }
    if(OPS_GetIntInput(&numData, &cDOF)) {
        opserr<<"WARNING invalid cDOF inputs\n";
        return -1;
    }
    cDOF--;
    double cc = 0.0;
    if (OPS_GetDouble(&numData, &cc) || cc == 0.0) {
        opserr<<"WARNING invalid ccoef inputs\n";
        return -1;
    }

    int rdf = numRemainingArgs / 3 - 1;
    ID rNode(rdf);
    ID rDOF(rdf);
    
     // constraint vector
    Vector Ccr(rdf);

    for(int i = 0; i < rdf; i++) {
        int rNodei, rDOFi;
        if(OPS_GetIntInput(&numData, &rNodei)) {
            opserr<<"WARNING invalid rNodeTag inputs\n";
            return -1;
        }
        if(OPS_GetIntInput(&numData, &rDOFi)) {
            opserr<<"WARNING invalid rDOF inputs\n";
            return -1;
        }
        double rci;
        if (OPS_GetDouble(&numData, &rci)) {
            opserr<<"WARNING invalid rcoef inputs\n";
            return -1;
        }
        rNode(i) = rNodei;
        rDOF(i) = rDOFi - 1;
        Ccr(i) = -rci / cc;
    }

    EQ_Constraint* theEQ = new EQ_Constraint(cNode,cDOF,Ccr,rNode,rDOF);
    if(theEQ == 0) {
	    opserr<<"WARNING: failed to create EQ_Constraint\n";
	    return -1;
    }
    if(theDomain->addEQ_Constraint(theEQ) == false) {
	    opserr<<"WARNING: failed to add EQ_Constraint to domain\n";
	    delete theEQ;
	    return -1;
    }
    return 0;
}

// constructor for FEM_ObjectBroker
EQ_Constraint::EQ_Constraint(int clasTag )		
:DomainComponent(nextTag++, clasTag),
 nodeRetained(0),nodeConstrained(0),constraint(0),constrDOF(0),retainDOF(0), initialized(false),
 dbTag1(0), dbTag2(0)
{
  numEQs++;
}

// constructor for Subclass
EQ_Constraint::EQ_Constraint(int nodeConstr, int constrainedDOF, 
			     ID &nodeRetain, ID &retainedDOF, int clasTag)
:DomainComponent(nextTag++, clasTag),
 nodeRetained(0), nodeConstrained(nodeConstr), 
 constraint(0), constrDOF(constrainedDOF), retainDOF(0), initialized(false), dbTag1(0), dbTag2(0)
{
    numEQs++;
  
    nodeRetained = new ID(nodeRetain);
    retainDOF = new ID(retainedDOF);    
    if (nodeRetained == 0 || nodeRetain.Size() != nodeRetained->Size() ||
            retainDOF == 0 || retainedDOF.Size() != retainDOF->Size()) { 
        opserr << "EQ_Constraint::EQ_Constraint - ran out of memory 1\n";
        exit(-1);
    }    

    // resize initial state
    Uc0 = 0.0;
    Ur0.resize(retainDOF->Size());
    Ur0.Zero();
}


// general constructor for ModelBuilder
EQ_Constraint::EQ_Constraint(int nodeConstr, int constrainedDOF, Vector &constr,
                                ID &nodeRetain, ID &retainedDOF)
:DomainComponent(nextTag++, CNSTRNT_TAG_EQ_Constraint), 
nodeRetained(0), nodeConstrained(nodeConstr), 
constraint(0), constrDOF(constrainedDOF), retainDOF(0), initialized(false), dbTag1(0), dbTag2(0)
{
    numEQs++;    
    nodeRetained = new ID(nodeRetain);
    retainDOF = new ID(retainedDOF);    
    if (nodeRetained == 0 || nodeRetain.Size() != nodeRetained->Size() ||
            retainDOF == 0 || retainedDOF.Size() != retainDOF->Size()) { 
        opserr << "EQ_Constraint::EQ_Constraint - ran out of memory 1\n";
        exit(-1);
    }    
    
    constraint = new Vector(constr);
    if (constraint == 0 || constraint->Size() != constr.Size()) { 
        opserr << "MP_Constraint::MP_Constraint - ran out of memory 2\n";
        exit(-1);
    }        
    
    // resize initial state
    Uc0 = 0.0;
    Ur0.resize(retainDOF->Size());
    Ur0.Zero();
}



EQ_Constraint::~EQ_Constraint()
{
    // invoke the destructor on the matrix and the two ID objects
    if (constraint != 0)
	delete constraint;
    if (nodeRetained != 0)
	delete nodeRetained;
    if (retainDOF != 0)
	delete retainDOF;    
    
    numEQs--;
    if (numEQs == 0)
        nextTag = 0;
}

void EQ_Constraint::setDomain(Domain* theDomain)
{
    // store initial state
    if (theDomain) {
        if (!initialized) {
            Node* theConstrainedNode = theDomain->getNode(nodeConstrained);
            if (theConstrainedNode == 0) {
                opserr << "FATAL EQ_Constraint::setDomain() - Constrained or Retained";
                opserr << " Node does not exist in Domain\n";
                opserr << nodeConstrained << endln;
                exit(-1);
            }
            const Vector& Uc = theConstrainedNode->getTrialDisp();
            int cdof = getConstrainedDOFs();
            if (cdof < 0 || cdof >= Uc.Size()) {
                opserr << "EQ_Constraint::setDomain FATAL Error: Constrained DOF " << cdof << " out of bounds [0-" << Uc.Size() - 1 << "]\n";
                exit(-1);
            }
            Uc0 = Uc(cdof);
            initialized = true;
            const ID& idr = getRetainedDOFs();
            for (int i = 0; i < nodeRetained->Size(); ++i) {
                Node* theRetainedNode = theDomain->getNode((*nodeRetained)(i));
                if (theRetainedNode == 0) {
                    opserr << "FATAL EQ_Constraint::setDomain() - Constrained or Retained";
                    opserr << " Node does not exist in Domain\n";
                    opserr << nodeRetained << endln;
                    exit(-1);
                }
                const Vector& Ur = theRetainedNode->getTrialDisp();
                int rdof = idr(i);
                if (rdof < 0 || rdof >= Ur.Size()) {
                    opserr << "EQ_Constraint::setDomain FATAL Error: Retained DOF " << rdof << " out of bounds [0-" << Ur.Size() - 1 << "]\n";
                    exit(-1);
                }
                Ur0(i) = Ur(rdof);
            }
        }
    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

const ID &
EQ_Constraint::getNodeRetained(void) const
{
    if (nodeRetained == 0) {
        opserr << "EQ_Constraint::getNodeRetained - no ID was set, ";
        opserr << "was recvSelf() ever called? or subclass incorrect?\n";	
        exit(-1);
    }

    // return the ID corresponding to retained nodes
    return *nodeRetained;    
}

int
EQ_Constraint::getNodeConstrained(void) const
{
    // return id of constrained node    
    return nodeConstrained;
}


int
EQ_Constraint::getConstrainedDOFs(void) const
{
    // return the ID corresponding to constrained DOF of Ccr
    return constrDOF;    
}


const ID &
EQ_Constraint::getRetainedDOFs(void) const
{
    if (retainDOF == 0) {
        opserr << "EQ_Constraint::getRetainedDOFs - no ID was set\n ";
        opserr << "was recvSelf() ever called? or subclass incorrect?\n";		
        exit(-1);
    }

    // return the ID corresponding to retained DOF of Ccr
    return *retainDOF;    
}

int 
EQ_Constraint::applyConstraint(double timeStamp)
{
    // does nothing EQ_Constraint objects are time invariant
    return 0;
}

bool
EQ_Constraint::isTimeVarying(void) const
{
    return false;
}


const Vector &
EQ_Constraint::getConstraint(void)
{
    if (constraint == 0) {
        opserr << "EQ_Constraint::getConstraint - no constraint was set\n";
        exit(-1);
    }    

    // return the constraint vector Ccr
    return *constraint;    
}

double EQ_Constraint::getConstrainedDOFsInitialDisplacement(void) const
{
    return Uc0;
}

const Vector& EQ_Constraint::getRetainedDOFsInitialDisplacement(void) const
{
    return Ur0;
}

int 
EQ_Constraint::sendSelf(int cTag, Channel &theChannel)
{
    static ID data(10);
    static Vector dataUc0(1);
    int dataTag = this->getDbTag();

    data(0) = this->getTag(); 
    data(1) = nodeConstrained;
    data(2) = constrDOF;
    if (constraint == 0) data(3) = 0; else data(3) = constraint->Size();
    if (nodeRetained== 0) data(4) = 0; else data(4) = nodeRetained->Size();    
    if (retainDOF == 0) data(5) = 0; else data(5) = retainDOF->Size();        
    
    // need two database tags for ID objects
    if (nodeRetained != 0 && dbTag1 == 0) 
      dbTag1 = theChannel.getDbTag();
    if (retainDOF != 0 && dbTag2 == 0) 
      dbTag2 = theChannel.getDbTag();

    data(6) = dbTag1;
    data(7) = dbTag2;
    data(8) = nextTag;
    data(9) = static_cast<int>(initialized);

    int result = theChannel.sendID(dataTag, cTag, data);
    if (result < 0) {
        opserr << "WARNING EQ_Constraint::sendSelf - error sending ID data\n";
        return result;
    }
    
    dataUc0(0) = Uc0;

    result = theChannel.sendVector(dataTag, cTag, dataUc0);
    if (result < 0) {
        opserr << "WARNING EQ_Constraint::sendSelf - error sending Vector data\n";
        return result;
    }

    if (constraint != 0 && constraint->Size() != 0) {
        int result = theChannel.sendVector(dataTag, cTag, *constraint);
        if (result < 0) {
            opserr << "WARNING EQ_Constraint::sendSelf ";
            opserr << "- error sending Vector data\n"; 
            return result;  
        }
    }

    if (nodeRetained != 0 && nodeRetained->Size() != 0) {
        int result = theChannel.sendID(dbTag1, cTag, *nodeRetained);
        if (result < 0) {
            opserr << "WARNING EQ_Constraint::sendSelf ";
            opserr << "- error sending nodeRetained data\n"; 
            return result;  
        }
    }

    if (retainDOF != 0 && retainDOF->Size() != 0) {
        int result = theChannel.sendID(dbTag2, cTag, *retainDOF);
        if (result < 0) {
            opserr << "WARNING EQ_Constraint::sendSelf ";
            opserr << "- error sending retainDOF data\n"; 
            return result;  
        }
    }
    
    // send initial displacement vectors at retained node.
    if (Ur0.Size() > 0) {
        int result = theChannel.sendVector(dbTag1, cTag, Ur0);
        if (result < 0) {
            opserr << "WARNING EQ_Constraint::sendSelf ";
            opserr << "- error sending retained initial displacement\n";
            return result;
        }
    }
    return 0;
}


int 
EQ_Constraint::recvSelf(int cTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    static ID data(10);
    static Vector dataUc0(1);
    int result = theChannel.recvID(dataTag, cTag, data);
    if (result < 0) {
        opserr << "WARNING EQ_Constraint::recvSelf - error receiving ID data\n";
        return result;  
    }    

    this->setTag(data(0));
    nodeConstrained = data(1);
    constrDOF = data(2);
    dbTag1 = data(6);
    dbTag2 = data(7);
    nextTag = data(8);
    initialized = static_cast<bool>(data(9));

    result = theChannel.recvVector(dataTag, cTag, dataUc0);
    if (result < 0) {
        opserr << "WARNING EQ_Constraint::recvSelf - error receiving Vector data\n";
        return result;  
    }

    Uc0 = dataUc0(0);

    int size = data(3); 
    if (size != 0) {
        constraint = new Vector(size);
        
        int result = theChannel.recvVector(dataTag, cTag, *constraint);
        if (result < 0) {
            opserr << "WARNING EQ_Constraint::recvSelf ";
            opserr << "- error receiving Vector data\n"; 
            return result;  
        }
    }    
    size = data(4);
    if (size != 0) {
        nodeRetained = new ID(size);
        int result = theChannel.recvID(dbTag1, cTag, *nodeRetained);
        if (result < 0) {
            opserr << "WARNING EQ_Constraint::recvSelf ";
            opserr << "- error receiving nodeRetained data\n"; 
            return result;  
        }	
    }
    
    size = data(5);
    if (size != 0) {
        retainDOF = new ID(size);
        int result = theChannel.recvID(dbTag2, cTag, *retainDOF);
        if (result < 0) {
            opserr << "WARNING EQ_Constraint::recvSelf ";
            opserr << "- error receiving retained data\n"; 
            return result;  
        }	
    }    
    
    // recv initial displacement vectors.
    if (retainDOF && retainDOF->Size() > 0)
        Ur0.resize(retainDOF->Size());
    else
        Ur0 = Vector();
    if (Ur0.Size() > 0) {
        int result = theChannel.recvVector(dbTag1, cTag, Ur0);
        if (result < 0) {
            opserr << "WARNING EQ_Constraint::recvSelf ";
            opserr << "- error receiving retained initial displacement\n";
            return result;
        }
    }
    return 0;
}



void
EQ_Constraint::Print(OPS_Stream &s, int flag)
{     
    s << "EQ_Constraint: " << this->getTag() << "\n";
    s << "Constrained Node: " << nodeConstrained << "at DOF: " << constrDOF + 1 << "\n";
    if (nodeRetained != 0 && retainDOF != 0)
        for (int i = 0; i < (*nodeRetained).Size(); i++)
            s << "Retained Node: " << (*nodeRetained)(i) << "at DOF: " << (*retainDOF)(i) + 1 << "\n";
    if (constraint != 0) {
        s << " Constraint vector: " << *constraint << "\n";
    }
}


