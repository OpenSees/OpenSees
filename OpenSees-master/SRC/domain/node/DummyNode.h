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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/node/DummyNode.h,v $
                                                                        
                                                                        
// File: ~/domain/node/DummyNode.h
//
// Written: fmk 
// Created: Fri Sep 20 15:34:47: 1996
// Revision: A
//
// Purpose: This file contains the class definition for DummyNode.
// DummyNodes are a type of node created and used by Subdomains for their
// exterior nodes. They reference a real node and most methods invoked on
// them are in turn invoked by the dummy node on the real node. The calls
// asking the real node to change its current state are ignored. The
// calls involving DOF\_Group are handled by the dummy node. \\
//
// 
// What: "@(#) DummyNode.h, revA"

#ifndef DummyNode_h
#define DummyNode_h

#include <Node.h>

class DOF_Group;

class DummyNode : public Node
{
  public:
    DummyNode();    
    DummyNode(Node &theRealNode);
    
    ~DummyNode();
    
    Node *getNode(void) const;
    
    void setDOF_GroupPtr(DOF_Group *theDOF_Grp);
    DOF_Group *getDOF_GroupPtr(void);
    int  getNumberDOF(void) const;

    const Matrix &getMass(void);
    int setMass(const Matrix &);
    
    const Vector &getCrds(void) const;

    const Vector &getDisp(void) ;
    const Vector &getVel(void) ;
    const Vector &getAccel(void) ;    

    const Vector &getTrialDisp(void) ;    
    const Vector &getTrialVel(void) ;    
    const Vector &getTrialAccel(void) ;    
    
    int setTrialDisp(const Vector &);    
    int setTrialVel(const Vector &);    
    int setTrialAccel(const Vector &);        

    int setIncrTrialDisp(const Vector &);    
    int setIncrTrialVel(const Vector &);    
    int setIncrTrialAccel(const Vector &);        
    
    void addUnbalancedLoad(const Vector &);        
    const Vector &getUnbalancedLoad(void) const;    
    void zeroUnbalancedLoad(void);
    int commitState();

    void Print(OPS_Stream &s) const;
    friend OPS_Stream &operator<<(OPS_Stream &s, const DummyNode &N);

    int addElementPtr(Element *);

    void setColor(int newColor);
    int getColor(void) const;

    int sendSelf(Channel &theChannel, FEM_ObjectBroker &theBroker);
    int recvSelf(Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
  private:
    Node *theRealNode;
    DOF_Group *theDOFGroup;
};

#endif

