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
// $Date: 2003-02-25 23:33:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthSection.h,v $
                                                                        
// Written: MHS
// Created: Sept 2000
//
// Description: This file contains the class definition for ZeroLengthSection.
// A ZeroLengthSection element is defined by two nodes with the same coordinate.
// A SectionForceDeformation object is associated with the nodes to
// provide the basic force-deformation relationship for the element.

#ifndef ZeroLengthSection_h
#define ZeroLengthSection_h

#include <Element.h>
#include <Matrix.h>

// Tolerance for zero length of element
#define	LENTOL 1.0e-6

class Node;
class Channel;
class SectionForceDeformation;
class Response;

class ZeroLengthSection : public Element
{
  public:
    
    ZeroLengthSection(int tag, 			      
	       int dimension,
	       int Nd1, int Nd2, 
	       const Vector& x,
	       const Vector& yprime,
		   SectionForceDeformation& theSection);

    ZeroLengthSection();    
    ~ZeroLengthSection();

    // public methods to obtain inforrmation about dof & connectivity    
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);	
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);    

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);    
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, Information &eleInformation);
    int getResponse(int responseID, Information &eleInformation);
    
  protected:
    
  private:
    // private methods
    void setUp (int Nd1, int Nd2, const Vector& x, const Vector& y);
    void setTransformation(void);
	void computeSectionDefs(void);

    // private attributes - a copy for each object of the class
    ID  connectedExternalNodes;         // contains the tags of the end nodes
    int dimension;                      // = 2 or 3 dimensions
    int numDOF;	                        // number of dof for ZeroLengthSection
    Matrix transformation;		// transformation matrix for orientation
	
    Matrix *A;	// Transformation matrix ... e = A*(u2-u1)
    Vector *v;	// Section deformation vector, the element basic deformations
    
    Matrix *K;	// Pointer to element stiffness matrix
    Vector *P;	// Pointer to element force vector
    
    Node *theNodes[2];
    
    SectionForceDeformation *theSection;	// Pointer to section object
    int order;		// Order of the section model
    
    // Class wide matrices for return
    static Matrix K6;
    static Matrix K12;
    
    // Class wide vectors for return
    static Vector P6;
    static Vector P12;
};

#endif




