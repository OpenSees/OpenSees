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
                                                                        
// $Revision$
// $Date$
// $Source$
                                                                        
                                                                        
#ifndef TPB1D_h
#define TPB1D_h

// Written: Troy/Fenz/Misra
// C++ Conversion by: fmk
// Created: 01/11
//
// Description: This file contains the class definition for TPB1D.
//
// What: "@(#) TPB1D.h, revA"

#include <Element.h>
#include <Matrix.h>

// Type of dimension of element NxDy has dimension x=1,2,3 and
// y=2,4,6,12 degrees-of-freedom for the element
enum Etype { D1N2, D2N4, D2N6, D3N6, D3N12 };

class Node;
class Channel;
class UniaxialMaterial;
class Response;

class TPB1D : public Element
{
  public:
    
  // Constructor for a single 1d material model
  TPB1D(int tag, 			      
	int Nd1, 
	int Nd2, 
	int dir,
	double *mu,
	double *R,
	double *h,
	double *D,
	double *d,
	double W);
  
  TPB1D();    
  ~TPB1D();
  
  const char *getClassType(void) const {return "TPB1D";};
  
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
  const Matrix &getDamp(void);
  const Matrix &getMass(void);
  
  void zeroLoad(void);	
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);    
  
  const Vector &getResistingForce(void);
  const Vector &getResistingForceIncInertia(void);            
  
  // public methods for element output
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);    
  int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
  
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);
  
  void updateDir (const Vector& x, const Vector& y);
  
 protected:
  
 private:
  // private attributes - a copy for each object of the class
  ID  connectedExternalNodes;         // contains the tags of the end nodes
  int dimension;                      // = 1, 2, or 3 dimensions
  int numDOF;	                      // number of dof for TPB1D
  int direction;
  
  Node *theNodes[2];

  double mu[3];
  double R[3];
  double h[3];
  double D[3];
  double d[3];
  double W;

  Matrix *theMatrix; 	    	// pointer to objects matrix (a class Matrix)
  Vector *theVector;      	// pointer to objects vector (a class Vector)
  
  UniaxialMaterial *theMaterial;      // array of pointers to 1d materials
  
  // vector pointers to initial disp and vel if present
  Vector *d0;
  
  // static data - single copy for all objects of the class	
  static Matrix TPB1DM2;   // class wide matrix for 2*2
  static Matrix TPB1DM4;   // class wide matrix for 4*4
  static Matrix TPB1DM6;   // class wide matrix for 6*6
  static Matrix TPB1DM12;  // class wide matrix for 12*12
  static Vector TPB1DV2;   // class wide Vector for size 2
  static Vector TPB1DV4;   // class wide Vector for size 4
  static Vector TPB1DV6;   // class wide Vector for size 6
  static Vector TPB1DV12;  // class wide Vector for size 12
};

#endif




