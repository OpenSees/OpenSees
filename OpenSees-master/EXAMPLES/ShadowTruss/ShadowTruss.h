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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-10-15 00:38:07 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/ShadowTruss/ShadowTruss.h,v $
                                                                        
#ifndef ShadowTruss_h
#define ShadowTruss_h

// Written: fmk 
// Created: 08/03
//
// Description: This file contains the interface for the ShadowTruss class.
//
// What: "@(#) ShadowTruss.h, revA"

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <Shadow.h>

class Node;
class Channel;
class UniaxialMaterial;
class ElementalLoad;

#define ELE_TAG_ShadowTruss 9001

class ShadowTruss : public Element, public Shadow
{
  public:
    // constructors
    ShadowTruss(int tag, 
		int Nd1, int Nd2, 
		UniaxialMaterial &theMaterial,
		double A, 
		double M, 
		Channel &theChannel,
		FEM_ObjectBroker &theObjectBroker);

    // destructor
    ~ShadowTruss();
    
    // public methods to obtain information about dof & connectivity
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);    

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
    const Vector &getResistingForce(void);

    // public methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);    
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, Information &eleInfo);
    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:
    ID msgData;

    // private member functions - only available to objects of the class
    double computeCurrentStrain(void) const;

    // private attributes - a copy for each object of the class
    UniaxialMaterial *theMaterial;       // pointer to a material
    ID  externalNodes;          	 // contains the id's of end nodes
    Matrix trans;       // hold the transformation matrix, could use a Vector
	                // if ever bother to change the Vector interface for
			// x-product.

    double L;		// length of truss (undeformed configuration) - set in setDomain()
    double A; 		// area of truss
    double M; 		// mass per unit volume of truss
    Node *end1Ptr, *end2Ptr; // two pointers to the nodes - these are set in setDomain()

    Vector *theLoad;     // pointer to the load vector P

    // static data - single copy for all objects of the class
    static Matrix trussK;   // class wide matrix for returning stiffness
    static Matrix trussM;   // class wide matrix for returning mass 
    static Vector trussR;   // class wide vector for returning residual
    static Node *theNodes[2];
};
#endif
