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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:00:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/node/DummyNode.cpp,v $
                                                                        
                                                                        
// File: ~/model/DummyNode.C
//
// Written: fmk 11/95
// Revised:
//
// This file contains the methods required for class DummyNode.
   
#include "DummyNode.h"
#include <classTags.h>

DummyNode::DummyNode()
:Node(0, NOD_TAG_DummyNode), 
 theRealNode(0), theDOFGroup(0)
{

}

DummyNode::DummyNode(Node &theNode)
:Node(theNode.getTag(), NOD_TAG_DummyNode), 
 theRealNode(&theNode), theDOFGroup(0)
{
    theRealNode = &theNode;
}



DummyNode::~DummyNode()
{

}

Node *
DummyNode::getNode() const
{

    return theRealNode;
}

void
DummyNode::setDOF_GroupPtr(DOF_Group *theDOF_Grp)
{
      theDOFGroup = theDOF_Grp;
}



DOF_Group *
DummyNode::getDOF_GroupPtr(void)
{
    return theDOFGroup;
}



int
DummyNode::getNumberDOF(void) const
{
    return theRealNode->getNumberDOF();
}





const Matrix &
DummyNode::getMass(void) 
{
    
    return theRealNode->getMass();
}




int
DummyNode::setMass(const Matrix &m)
{
    return theRealNode->setMass(m);
}



const Vector &
DummyNode::getCrds() const
{
    return theRealNode->getCrds();
}


const Vector &
DummyNode::getDisp(void) 
{
    return theRealNode->getDisp();
}


const Vector &
DummyNode::getVel(void) 
{
    return theRealNode->getVel();
}


const Vector &
DummyNode::getAccel(void) 
{
    return theRealNode->getAccel();
}


const Vector &
DummyNode::getTrialDisp(void) 
{
    return theRealNode->getTrialDisp();
}



const Vector &
DummyNode::getTrialVel(void) 
{
    return theRealNode->getTrialVel();
}



const Vector &
DummyNode::getTrialAccel(void) 
{
    return theRealNode->getTrialAccel();
}





int
DummyNode::setTrialDisp(const Vector &displ)
{
    return 0;
}

int
DummyNode::setTrialVel(const Vector &vel)
{
    return 0;
}


int
DummyNode::setTrialAccel(const Vector &accel)
{
    return 0;
}


int
DummyNode::setIncrTrialDisp(const Vector &incrDispl)
{
    return 0;
}


int
DummyNode::setIncrTrialVel(const Vector &incrVel)
{
    return 0;
}


int
DummyNode::setIncrTrialAccel(const Vector &incrAccel)
{
    return 0;
}


int
DummyNode::commitState()
{
    return 0;
}




void 
DummyNode::addUnbalancedLoad(const Vector &add)
{
    theRealNode->addUnbalancedLoad(add);
}



const Vector &
DummyNode::getUnbalancedLoad(void) const
{
    return theRealNode->getUnbalancedLoad();
}


void 
DummyNode::zeroUnbalancedLoad(void)
{
    theRealNode->zeroUnbalancedLoad();
}



int
DummyNode::addElementPtr(Element *elePtr)
{
//    return theRealNode->addElementPtr(elePtr);
  return 0;
}



void
DummyNode::setColor(int newColor)
{
//    theRealNode->setColor(newColor);
}

int
DummyNode::getColor(void) const
{
//    return theRealNode->getColor();
  return 0;
}

int 
DummyNode::sendSelf(Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    opserr << "DummyNode::sendSelf - should never be called\n";
    opserr << "sending acual node\n";    
    
    return theRealNode->sendSelf(theChannel,theBroker);
}

int 
DummyNode::recvSelf(Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    opserr << "DummyNode::recvSelf - should never be called\n";
    return 0;
}    

void
DummyNode::Print(OPS_Stream &s) const
{
    theRealNode->Print(s);
}
  
OPS_Stream &operator<<(OPS_Stream &s, const DummyNode &N)
{
    N.Print(s);
    return s;
}


