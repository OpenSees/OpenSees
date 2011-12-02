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
// $Date: 2005-10-20 23:52:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/beam2d/beam2d03.h,v $
                                                                        
                                                                        
// File: ~/element/beam2d03.h
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for beam2d03.
// beam2d03 is a plane frame member.

#ifndef beam2d03_h
#define beam2d03_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;

class beam2d03 : public Element
{
  public:
    beam2d03();        
    beam2d03(int tag, double A, double E, double I, int Nd1, int Nd2);
    ~beam2d03();

    void setDomain(Domain *theDomain);
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    int revertToLastCommit(void);        
    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag =0);    
    
  private:
    const Matrix &getStiff(void);    
    double A,E,I;
    double L,sn,cs;
    ID  connectedExternalNodes;    
    Node *theNodes[2];

    Matrix k; // the stiffness matrix
	
    Vector rForce;
    Vector load;
    Matrix trans;  // hold part of transformation matrix
};

#endif


