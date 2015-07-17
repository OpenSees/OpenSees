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
// $Date: 2013-05-31 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/multipleShearSpring/MultipleShearSpring.h,v $

// Written: Ken Ishii
// Created: June 2012
//
// Multiple Shear Spring (MSS) model
//
// Description: This file contains the class definition for MultipleShearSpring.

#ifndef MultipleShearSpring_h
#define MultipleShearSpring_h

#include <Element.h>
#include <Matrix.h>

class Channel;
class UniaxialMaterial;
class Response;

class MultipleShearSpring : public Element
{
 public:
  // constructor
  MultipleShearSpring(int Tag, int Nd1, int Nd2,
		      int NSpring,
		      UniaxialMaterial *Material,
		      double LimDisp,
		      const Vector OriYp, const Vector OriX = 0,
		      double Mass = 0.0);

  MultipleShearSpring(int Tag, int Nd1, int Nd2,
		      UniaxialMaterial **theMaterials,
		      int NSpring,
		      double LimDisp,
		      const Vector OriYp, const Vector OriX = 0,
		      double Mass = 0.0);

  MultipleShearSpring();
  
  // destructor
  ~MultipleShearSpring();
  
  // method to get class type
  const char *getClassType() const {return "MultipleShearSpring";};
  
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
  const Matrix &getMass();
  
  void zeroLoad();
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);
  
  const Vector &getResistingForce();
  const Vector &getResistingForceIncInertia();
  
  // public methods for element output
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);

  void Print(OPS_Stream &s, int flag = 0);    
  
  // public methods for element recorder
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInfo);
  
 protected:
  
 private:
  // private methods
  void setUp();
  
  // private attributes - a copy for each object of the class
  ID connectedExternalNodes;        // contains the tags of the end nodes
  Node *theNodes[2];                // array of nodes
  UniaxialMaterial **theMaterials;  // materials
  
  // parameters
  int nSpring; //number of shear springs in MSS
  double *cosTht; //arrangement of each spring (cos)
  double *sinTht; //arrangement of each spring (sin)
  Vector oriX;   // local x direction
  Vector oriYp;  // local yp direction
  double mass; // mass of element

  // calculation of Feq and Seq
  double limDisp; //minimum deformation to calculate Feq and Seq (if limDisp is 0, never calculate)
  UniaxialMaterial *dmyMssMaterial; //imaginary material to calculate Feq and Seq
  double mssFeq; //equivalent coefficient for force
  double mssSeq; //equivalent coefficient for stiffness

  
  // transformation
  Matrix Tgl; // transformation matrix from global to local system
  Matrix Tlb; // transformation matrix from local to basic system
  
  // displacement, force, stiffness
  Vector basicDisp;  //ub
  Vector localDisp;  //ul
  Vector basicForce; //qb
  Matrix basicStiff; //kb  
  Matrix basicStiffInit; //kbInit

  static Matrix theMatrix;
  static Vector theVector;
  static Vector theLoad;
};

#endif
