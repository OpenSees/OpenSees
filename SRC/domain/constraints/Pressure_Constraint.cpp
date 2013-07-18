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

Pressure_Constraint::Pressure_Constraint(int classTag)
    :DomainComponent(0,classTag), pTag(0), fluidEleTags(), otherEleTags(), gravity(0)
{
}

Pressure_Constraint::Pressure_Constraint(int classTag, int nodeId, int ptag)
    :DomainComponent(nodeId,classTag), pTag(ptag), fluidEleTags(), otherEleTags(), gravity(0)
{
}

Pressure_Constraint::Pressure_Constraint(int nodeId, int ptag)
    :DomainComponent(nodeId,CNSTRNT_TAG_Pressure_Constraint), 
     pTag(ptag), fluidEleTags(), otherEleTags(), gravity(0)
{
}

Pressure_Constraint::Pressure_Constraint(int classTag, int nodeId, double g)
    :DomainComponent(nodeId,classTag), pTag(nodeId), fluidEleTags(), otherEleTags(), gravity(g)
{
}

Pressure_Constraint::Pressure_Constraint(int nodeId, double g)
    :DomainComponent(nodeId,CNSTRNT_TAG_Pressure_Constraint), 
     pTag(nodeId), fluidEleTags(), otherEleTags(), gravity(g)
{
}

Pressure_Constraint::~Pressure_Constraint()
{
    Domain* theDomain = this->getDomain();
    if(theDomain != 0) {
        Node* pNode = theDomain->removeNode(pTag);
        if(pNode != 0) {
            delete pNode;
        }
    }
}

void
Pressure_Constraint::setDomain(Domain* theDomain)
{
    // get old domain
    // Domain* oldDomain = this->getDomain();
    // if(oldDomain != 0) {
    //     Node* pNode = oldDomain->removeNode(pTag);
    //     if(pNode != 0) {
    //         delete pNode;
    //     }
    // }

    // set new domain
    DomainComponent::setDomain(theDomain);

    // check domain
    if(theDomain == 0) return;
    int nodeId = this->getTag();
    Node* pNode = 0;
    if(pTag != nodeId) {
        pNode = theDomain->getNode(pTag);
        if(pNode != 0) return;
    } else {
        pTag = findNodeTag(theDomain);
    }

    // check the node
    Node* theNode = theDomain->getNode(nodeId);
    if(theNode == 0) {
        opserr<<"WARNING: node "<<nodeId<<" does not exist ";
        opserr<<"-- Pressure_Constraint::setDomain\n";
        return;
    }
    
    // check ndm
    const Vector& coord = theNode->getCrds();
    int ndm = coord.Size();
    if(ndm<2 || ndm>3) {
        opserr<<"WARNING: bad ndm "<<ndm<<" of node "<<nodeId;
        opserr<<" -- Pressure_Constraint::setDomain\n";
        return;
    }

    // check ndf
    int ndf = theNode->getNumberDOF();
    if(ndf < ndm) {
        opserr<<"WARNING: wrong ndf of node "<<nodeId;
        opserr<<": must be at least "<<ndm;
        opserr<<" -- Pressure_Constraint::setDomain\n";
        return;
    }

    // create pressure node
    if(ndm == 2) {
        pNode = new Node(pTag, ndm+1, coord(0), coord(1));
    } else {
        pNode = new Node(pTag, ndm+1, coord(0), coord(1), coord(2));
    }
    if(pNode == 0) {
        opserr<<"WARNING: run out of memory to create Node";
        opserr<<" -- Pressure_Constraint::setDomain\n";
        return;
    }
    if(theDomain->addNode(pNode) == false) {
        opserr<<"WARNING: failed to add node to domain";
        opserr<<" -- Pressure_Constraint::setDomain\n";
        delete pNode;
        return;
    }



}

int
Pressure_Constraint::getPressureNode()
{
    Domain* theDomain = this->getDomain();
    if(theDomain == 0) {
        opserr<<"WARNING: domain has not been set";
        opserr<<" -- Pressure_Constraint::getPressureNode\n";
        return 0;
    }
    return pTag;
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

void 
Pressure_Constraint::newStep(double dt, Vector& U, Vector& Udot, Vector& Udotdot)
{
    if(this->isIsolated()) {

        // check domain
        Domain* theDomain = this->getDomain();
        if(theDomain == 0) {
            opserr<<"WARNING: domain has not been set";
            opserr<<" -- Pressure_Constraint::applyConstraint\n";
            return;
        }

        // check the node
        int nodeId = this->getTag();
        Node* theNode = theDomain->getNode(nodeId);
        if(theNode == 0) {
            opserr<<"WARNING: node "<<nodeId<<" does not exist ";
            opserr<<"-- Pressure_Constraint::applyConstraint\n";
            return;
        }

        // DOF_group
        DOF_Group* theDOF = theNode->getDOF_GroupPtr();
        if(theDOF == 0) {
            opserr<<"WARNING: no DOF_Group is found with node "<<nodeId;
            opserr<<" -- Pressure_Constraint::newStep\n";
            return;
        }
        const ID& id = theDOF->getID();

        // get commiated states
        const Vector& vel = theNode->getVel();
        const Vector& disp = theNode->getDisp();
        const Vector& coord = theNode->getCrds();
        int ndm = coord.Size();

        // set new states to the ndm-1 direction
        if(id(ndm-1) >= 0) {
            Udotdot(id(ndm-1)) = gravity;
            Udot(id(ndm-1)) = vel(ndm-1) + gravity*dt;
            U(id(ndm-1)) = disp(ndm-1) + vel(ndm-1)*dt + 0.5*gravity*dt*dt;
        }
    }
}


int 
Pressure_Constraint::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;
    int dataTag = this->getDbTag();

    int sf = fluidEleTags.Size();
    int so = otherEleTags.Size();

    // send size of IDs
    static ID idsize(2);
    idsize(0) = sf;
    idsize(1) = so;
    res = theChannel.sendID(dataTag, commitTag, idsize);
    if(res < 0) {
        opserr<<"WARNING: Pressure_Constraint::sendSelf - "<<this->getTag()<<" failed to send idsize\n";
        return -1;
    }

    // send vector
    static Vector data(3+sf+so);
    data(0) = this->getTag();
    data(1) = pTag;
    data(2) = gravity;
    for(int i=0; i<sf; i++) {
        data(3+i) = fluidEleTags(i);
    }
    for(int i=0; i<so; i++) {
        data(3+sf+i) = otherEleTags(i);
    }
    res = theChannel.sendVector(dataTag, commitTag, data);
    if(res < 0) {
        opserr<<"WARNING: Pressure_Constraint::sendSelf - "<<this->getTag()<<" failed to send vector\n";
        return -2;
    }

    return 0;
}

int 
Pressure_Constraint::recvSelf(int commitTag, Channel &theChannel, 
                                  FEM_ObjectBroker &theBroker)
{
    int res = 0;
    int dataTag = this->getDbTag();

    // receive idsize
    static ID idsize(2);
    res = theChannel.recvID(dataTag, commitTag, idsize);
    if(res < 0) {
        opserr<<"WARNING: PFEMElement2D::recvSelf - failed to receive idsize\n";
        return -1;
    }
    int sf = idsize(0);
    int so = idsize(1);
    fluidEleTags.resize(sf);
    otherEleTags.resize(so);

    // receive vector
    static Vector data(3+sf+so);
    res = theChannel.recvVector(dataTag, commitTag, data);
    if(res < 0) {
        opserr<<"WARNING: Pressure_Constraint::recvSelf - "<<this->getTag()<<" failed to receive vector\n";
        return -2;
    }
    this->setTag((int)data(0));
    pTag = (int)data(1);
    gravity = data(2);
    for(int i=0; i<sf; i++) {
        fluidEleTags(i) = (int)data(3+i);
    }
    for(int i=0; i<so; i++) {
        otherEleTags(i) = (int)data(3+sf+i);
    }
    return 0;
}

void 
Pressure_Constraint::Print(OPS_Stream& s, int flag)
{
    s << "Pressure_Constraint: " << this->getTag() << "\n";
    s << "pressure node -- "<<pTag<<"\n";
}

int 
Pressure_Constraint::findNodeTag(Domain* theDomain)
{
    NodeIter& theNodes = theDomain->getNodes();
    Node* firstNode = theNodes();
    if(firstNode == 0) {
        return -1;
    }
    int tag = firstNode->getTag();
    return tag-1;
}

void
Pressure_Constraint::setGravity(double g)
{
    gravity = g;
}
