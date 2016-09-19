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
    :DomainComponent(0,classTag), pTag(0), fluidEleTags(), otherEleTags()
{
}

Pressure_Constraint::Pressure_Constraint(int classTag, int nodeId, int ptag, int ndf)
    :DomainComponent(nodeId,classTag), pTag(ptag), fluidEleTags(), otherEleTags(), pndf(ndf)
{
}

Pressure_Constraint::Pressure_Constraint(int nodeId, int ptag, int ndf)
    :DomainComponent(nodeId,CNSTRNT_TAG_Pressure_Constraint), 
     pTag(ptag), fluidEleTags(), otherEleTags(), pndf(ndf)
{
}

Pressure_Constraint::Pressure_Constraint(int nodeId, int ndf)
    :DomainComponent(nodeId,CNSTRNT_TAG_Pressure_Constraint), 
     pTag(nodeId), fluidEleTags(), otherEleTags(), pndf(ndf)
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

    // check pndf
    int nodeId = this->getTag();
    if(pndf <= 0) {
        return;
    }

    // create pressure node
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
        pNode = new Node(pTag, pndf, coord(0), coord(1));
    } else {
        pNode = new Node(pTag, pndf, coord(0), coord(1), coord(2));
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

Node*
Pressure_Constraint::getPressureNode()
{
    if(pTag == this->getTag()) return 0;
    if(pndf <= 0) return 0;
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
    Node* pnode = this->getPressureNode();
    if (pnode == 0) return;
    const Vector& vel = pnode->getVel();
    Vector newvel(vel);
    newvel.Zero();
    newvel(0) = p;
    pnode->setTrialVel(newvel);
    pnode->commitState();
}

double
Pressure_Constraint::getPressure(int last)
{
    if(pndf <= 0) return 0.0;
    Domain* theDomain = this->getDomain();
    if(theDomain == 0) {
        opserr<<"WARNING: domain has not been set";
        opserr<<" -- Pressure_Constraint::getPressureNode\n";
        return 0;
    }
    int ntag = this->getTag();
    if(ntag == pTag) return 0.0;
    Node* pNode = theDomain->getNode(pTag);
    if(pNode == 0) return 0.0;
    const Vector& vel = pNode->getVel();
    if(last == 1) {
        if(vel.Size()==0) return 0.0;
        return vel(0);
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
    int res = 0;
    int dataTag = this->getDbTag();

    // send size of IDs
    ID idsize(3);
    idsize(0) = pTag;
    idsize(1) = pndf;
    idsize(2) = this->getTag();
    res = theChannel.sendID(dataTag, commitTag, idsize);
    if(res < 0) {
        opserr<<"WARNING: Pressure_Constraint::sendSelf - ";
        opserr<<this->getTag()<<" failed to send idsize\n";
        return -1;
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
    ID idsize(3);
    res = theChannel.recvID(dataTag, commitTag, idsize);
    if(res < 0) {
        opserr<<"WARNING: PFEMElement2D::recvSelf - failed to receive idsize\n";
        return -1;
    }
    pTag = idsize(0);
    pndf = idsize(1);
    this->setTag(idsize(2));

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

