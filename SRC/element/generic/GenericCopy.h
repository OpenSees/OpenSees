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

#ifndef GenericCopy_h
#define GenericCopy_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 11/06
// Revision: A
//
// Description: This file contains the class definition for GenericCopy.
// GenericCopy is an exact copy of an already existing element in the
// domain. As such it returns the same stiffnesses and resisting forces
// as the source element. This is useful if the source element is a
// genericClient element that communicates with an experimental element.

#include <Element.h>
#include <Matrix.h>


class GenericCopy : public Element
{
public:
    // constructors
    GenericCopy(int tag, ID nodes, int srcTag);
    GenericCopy();
    
    // destructor
    ~GenericCopy();
    
    // method to get class type
    const char *getClassType() const {return "GenericCopy";};
    
    // public methods to obtain information about dof & connectivity
    int getNumExternalNodes() const;
    const ID &getExternalNodes();
    Node **getNodePtrs();
    int getNumDOF();
    void setDomain(Domain *theDomain);
    
    // public methods to set the state of the element
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();
    
    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getDamp();
    const Matrix &getMass();
    
    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
    
    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();
    
    // public methods for element output
    int sendSelf(int commitTag, Channel &sChannel);
    int recvSelf(int commitTag, Channel &rChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
    void Print(OPS_Stream &s, int flag = 0);
    
    // public methods for element recorder
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);
    
protected:
    
private:
    // private attributes - a copy for each object of the class
    ID connectedExternalNodes;  // contains the tags of the end nodes
    
    int numExternalNodes;       // number of external nodes
    int numDOF;                 // number of total DOF
    
    int srcTag;                 // tag of the source element
    Element *theSource;         // pointer to the source element
    
    Matrix theMatrix;           // objects matrix
    Vector theVector;           // objects vector
    Vector theLoad;             // load vector
    Matrix theInitStiff;        // initial stiffness matrix
    Matrix theMass;             // mass matrix

    bool initStiffFlag;
    bool massFlag;
    
    Node **theNodes;
};

#endif
