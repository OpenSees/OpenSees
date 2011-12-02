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
                                                                        
// $Revision: 1.8 $
// $Date: 2006-08-04 19:13:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/truss/CorotTruss.h,v $

#ifndef CorotTruss_h
#define CorotTruss_h

// Written: MHS 
// Created: May 2001
//
// Description: This file contains the class definition for CorotTruss,
// a small strain, large displacement corotational space truss element,
// as described by Crisfield in "Nonlinear Finite Element Analysis of
// Solids and Structures", Vol. 1, 1991, J.T. Wiley.

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Node;
class Channel;
class UniaxialMaterial;

class CorotTruss : public Element
{
  public:
    CorotTruss(int tag, int dim,
	  int Nd1, int Nd2, 
	  UniaxialMaterial &theMaterial,
	  double A, double rho=0.0);
    
    CorotTruss();    
    ~CorotTruss();

    const char *getClassType(void) const {return "CorotTruss";};

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
    int update(void);
    
    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);    

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

    Response *setResponse(const char **argv, int argc, Information &eleInfo, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:
   
    // private attributes - a copy for each object of the class
    UniaxialMaterial *theMaterial;  // pointer to a material
    ID  connectedExternalNodes;     // contains the tags of the end nodes
    int numDOF;	                    // number of dof for CorotTruss
    int numDIM;                     // number of dimensions

    double Lo;	        // initial length of truss
    double Ln;		// current length of truss
    double d21[3];	// current displacement offsets in basic system
    double A; 	        // area of CorotTruss
    double rho; 	        // mass density per unit length

    Node *theNodes[2];

    Matrix R;	// Rotation matrix

    Matrix *theMatrix;

    static Matrix M2;
    static Matrix M4;
    static Matrix M6;
    static Matrix M12;
    
    Vector *theVector;

    static Vector V2;
    static Vector V4;
    static Vector V6;
    static Vector V12;
};

#endif




