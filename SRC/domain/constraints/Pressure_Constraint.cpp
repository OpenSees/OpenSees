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
// $Date: 2012-08-21 12:59:06 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/Pressure_Constraint.h,v $

// Written: Minjie
// Created: 08/12
// Revision: A                                                             
//
// Purpose: This file contains the implementation of class Pressure_Constraint.
//
// The class Pressure_Constraint interface:
//

#include <Pressure_Constraint.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Node.h>
#include <NodeIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <Domain.h>
#include <Matrix.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <NodalLoad.h>
#include <Element.h>
#include <DOF_Group.h>
#include <elementAPI.h>

int
OPS_ADD_RUNTIME_IXV(OPS_Pressure_Constraint)
{
    Domain* theDomain = OPS_GetDomain();

    if (theDomain == 0) {
        opserr << "WARNING: domain is not defined\n";
        return -1;
    }

    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING: need nodeTag, pNodeTag\n";
        return -1;
    }

    // get node tags
    int tags[2];
    int num = 2;
    if (OPS_GetIntInput(&num, &tags[0]) < 0) {
        opserr << "WARNING: invalid node tag\n";
        return -1;
    }

    // set domain
    Pressure_Constraint* pc = new Pressure_Constraint(tags[0], tags[1]);
    if (pc == 0) {
        opserr << "WARNING: failed to create pc\n";
        return -1;
    }

    if (theDomain->addPressure_Constraint(pc) == false) {
        opserr << "WARNING: failed to add pc to domain\n";
        delete pc;
        return -1;
    }
    
    return 0;
}

Pressure_Constraint::Pressure_Constraint(int classTag)
    :DomainComponent(0,classTag), pTag(0), fluidEleTags(), otherEleTags(),
     pval(0), freesurf(false)
{
}

Pressure_Constraint::Pressure_Constraint(int nodeId, int ptag)
    :DomainComponent(nodeId,CNSTRNT_TAG_Pressure_Constraint),
     pTag(ptag), fluidEleTags(), otherEleTags(), pval(0), freesurf(false)
{
}

Pressure_Constraint::Pressure_Constraint(int nodeId, double val)
    :DomainComponent(nodeId,CNSTRNT_TAG_Pressure_Constraint),
     pTag(nodeId), fluidEleTags(), otherEleTags(),
     pval(0), freesurf(false)
{
    pval = new double[2];
    pval[0] = val;
    pval[1] = 0.0;
}

Pressure_Constraint::~Pressure_Constraint()
{
    Domain* theDomain = this->getDomain();
    if(theDomain != 0) {
	if (pval == 0) {
	    Node* pNode = theDomain->removeNode(pTag);
	    if(pNode != 0) {
		delete pNode;
	    }
	}
    }
    if (pval != 0) delete [] pval;
}

void
Pressure_Constraint::setDomain(Domain* theDomain)
{
    freesurf = false;
    // set new domain
    DomainComponent::setDomain(theDomain);

    // check domain
    if(theDomain == 0) return;

    // check the node
    int nodeId = this->getTag();
    Node* theNode = theDomain->getNode(nodeId);
    if(theNode == 0) {
        opserr<<"WARNING: node "<<nodeId<<" does not exist ";
        opserr<<"-- Pressure_Constraint::setDomain\n";
        return;
    }

    // check pval
    if (pval != 0) return;

    if (pTag == nodeId) {
	opserr << "WARNING: pressure node has the same tag as the PC\n";
	return;
    }

    // check the pressure node
    Node* pNode = theDomain->getNode(pTag);
    if(pNode == 0) {
        opserr<<"WARNING: pressure node "<<pTag<<" does not exist ";
        opserr<<"-- Pressure_Constraint::setDomain\n";
        return;
    }
}

Node*
Pressure_Constraint::getPressureNode()
{
    if(pval != 0) return 0;

    Domain* theDomain = this->getDomain();
    if(theDomain == 0) {
        opserr<<"WARNING: domain has not been set";
        opserr<<" -- Pressure_Constraint::getPressureNode\n";
        return 0;
    }
    return theDomain->getNode(pTag);
}

void
Pressure_Constraint::setPressure(double p)
{
    if (pval != 0) {
	pval[0] = p;
	return;
    }

    Node* pnode = this->getPressureNode();
    if (pnode == 0) return;
    const Vector& vel = pnode->getVel();
    Vector newvel(vel);
    newvel.Zero();
    newvel(0) = p;
    pnode->setTrialVel(newvel);
    pnode->commitState();
}

void
Pressure_Constraint::setPdot(double pdot)
{
    if (pval != 0) {
    	pval[1] = pdot;
    	return;
    }

    Node* pnode = this->getPressureNode();
    if (pnode == 0) return;
    const Vector& accel = pnode->getAccel();
    Vector newaccel(accel);
    newaccel.Zero();
    newaccel(0) = pdot;
    pnode->setTrialAccel(newaccel);
    pnode->commitState();
}

double
Pressure_Constraint::getPressure(int last)
{
    if (pval != 0) {
	return pval[0];
    }

    Domain* theDomain = this->getDomain();
    if(theDomain == 0) {
        opserr<<"WARNING: domain has not been set";
        opserr<<" -- Pressure_Constraint::getPressureNode\n";
        return 0;
    }
    Node* pNode = theDomain->getNode(pTag);
    if(pNode == 0) return 0.0;
    const Vector& vel = pNode->getVel();
    if(last == 1) {
        if(vel.Size()==0) return 0.0;
        return vel(0);
    }
    return 0.0;
}

double
Pressure_Constraint::getPdot(int last)
{
    if (pval != 0) {
	return pval[1];
    }

    Domain* theDomain = this->getDomain();
    if(theDomain == 0) {
        opserr<<"WARNING: domain has not been set";
        opserr<<" -- Pressure_Constraint::getPressureNode\n";
        return 0;
    }

    Node* pNode = theDomain->getNode(pTag);
    if(pNode == 0) return 0.0;
    const Vector& accel = pNode->getAccel();
    if(last == 1) {
        if(accel.Size()==0) return 0.0;
        return accel(0);
    }
    return 0.0;
}

const ID&
Pressure_Constraint::getFluidElements()
{
    Domain* theDomain = this->getDomain();
    if(theDomain == 0) {
        opserr<<"WARNING: domain has not been set";
        opserr<<" -- Pressure_Constraint::getConnectedElement\n";
    }
    return fluidEleTags;
}

const ID&
Pressure_Constraint::getOtherElements()
{
    Domain* theDomain = this->getDomain();
    if(theDomain == 0) {
        opserr<<"WARNING: domain has not been set";
        opserr<<" -- Pressure_Constraint::getConnectedElement\n";
    }
    return otherEleTags;
}

void
Pressure_Constraint::connect(int eleId, bool fluid)
{
    Domain* theDomain = this->getDomain();
    if(theDomain == 0) {
        opserr<<"WARNING: domain has not been set";
        opserr<<" -- Pressure_Constraint::connect\n";
        return;
    }

    Element* theEle = theDomain->getElement(eleId);
    if(theEle == 0) {
        opserr<<"WARNING: element "<<eleId<<" does not exist ";
        opserr<<"-- Pressure_Constraint::connect\n";
        return;
    }

    if(fluid) {
        fluidEleTags.insert(eleId);
    } else {
        bool isfluid = false;
        for(int i=0; i<fluidEleTags.Size(); i++) {
            if(fluidEleTags(i) == eleId) {
                isfluid = true;
                break;
            }
        }
        if(!isfluid) {
            otherEleTags.insert(eleId);
        }
    }
}

void
Pressure_Constraint::disconnect(int eleId)
{
    fluidEleTags.removeValue(eleId);
}

void
Pressure_Constraint::disconnect()
{
    otherEleTags = ID();
}

bool
Pressure_Constraint::isFluid() const
{
    return fluidEleTags.Size()>0 && otherEleTags.Size()==0;
}

bool
Pressure_Constraint::isInterface() const
{
    return fluidEleTags.Size()>0 && otherEleTags.Size()>0;
}

bool
Pressure_Constraint::isStructure() const
{
    return fluidEleTags.Size()==0 && otherEleTags.Size()>0;
}

bool
Pressure_Constraint::isIsolated() const
{
    return fluidEleTags.Size()==0 && otherEleTags.Size()==0;
}

int 
Pressure_Constraint::sendSelf(int commitTag, Channel& theChannel)
{
    return 0;
}

int 
Pressure_Constraint::recvSelf(int commitTag, Channel &theChannel, 
                                  FEM_ObjectBroker &theBroker)
{
    return 0;
}

void 
Pressure_Constraint::Print(OPS_Stream& s, int flag)
{
    s << "Pressure_Constraint: " << this->getTag() << "\n";
    s << "pressure node -- "<<pTag<<"\n";
}
