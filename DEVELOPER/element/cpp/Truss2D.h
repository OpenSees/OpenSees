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
// $Date: 2008/12/10 00:05:21 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewElement/cpp/Truss2D.h,v $
                                                                        
#ifndef Truss2D_h
#define Truss2D_h

// Written: fmk 
//
// Description: This file contains the interface for the Truss2D class.
// It defines the class interface and the class attributes. Truss2D 
// provides the abstraction of a simple truss element for 2-d problems.
// The stress-strain relationship for the truss being performed in the
// associated UniaxialMaterial object.
//
// What: "@(#) Truss2D.h, revA"

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class UniaxialMaterial;

class Truss2D : public Element
{
  public:
    // constructors
    Truss2D(int tag, 
	    int Nd1, int Nd2, 
	    UniaxialMaterial &theMaterial,
	    double A); 

    Truss2D();    
    
    // destructor
    ~Truss2D();

    
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

    // public methods to obtain stiffness
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);

    // public method to obtain resisting force
    const Vector &getResistingForce(void);

    // public methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:
    // private member functions - only available to objects of the class
    double computeCurrentStrain(void) const;
    
    // private attributes - a copy for each object of the class
    UniaxialMaterial *theMaterial;       // pointer to a material
    ID  externalNodes;          	 // contains the id's of end nodes
    Matrix trans;       // hold the transformation matrix
    double L;		// length of truss (undeformed configuration) - set in setDomain()
    double A; 		// area of truss
    Node *theNodes[2];  // node pointers

    // static data - single copy for all objects of the class
    static Matrix trussK;   // class wide matrix for returning stiffness
    static Vector trussR;   // class wide vector for returning residual
};
#endif

