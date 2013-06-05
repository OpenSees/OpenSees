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
// $Date: 2013-06-05 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/MSSWithMNS/MSSWithMNS.h,v $

// Written: Ken Ishii
// Created: Jan 2013
//
// MSSWithMNS model
//
// Description: This file contains the class definition for MSSWithMNS.
//

#ifndef MSSWithMNS_h
#define MSSWithMNS_h

#include <Element.h>
#include <Matrix.h>

class Channel;
class UniaxialMaterial;
class Response;

class MSSWithMNS : public Element
{
 public:
  // constructor
  MSSWithMNS(int Tag, int Nd1, int Nd2,
	     int Shape, double Size, double TotalRubber,
	     int NMSS, UniaxialMaterial *MatMSS, double LimDisp,
	     int NMNS, UniaxialMaterial *MatMNS, double Lambda,
	     const Vector OriYp, const Vector OriX = 0,
	     double Mass = 0.0);
  MSSWithMNS();
  
  // destructor
  ~MSSWithMNS();
  
  // method to get class type
  const char *getClassType() const {return "MSSWithMNS";};
  
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
  int displaySelf(Renderer &theViewer, int displayMode, float fact);    
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

  UniaxialMaterial **theINodeMNSMaterials;  // material‚ðnMNS^2ŒÂ iNode‘¤
  
  //shape
  int shape; // 1:circular, 2:square
  double size; //diameter (circular shape), length of edge (square shape)
  double totalRubber; //total rubber thickness
  double totalHeight; //total height of bearing
  //double hgt; //imaginary length of axial spring


  //MSS
  int nMSS; //number of shear springs in MSS
  double limDisp; //minimum deformation to calculate Feq and Seq (if limDisp is 0, never calculate)

  double *cosTht; //arrangement of each spring (cos)
  double *sinTht; //arrangement of each spring (sin)


  //MNS
  int nMNS; //section is divided into (nDivide)*(nDivide) springs
  double lambda; //parameter, =(D/t)*sqrt(3G/K)
  double incA; //area of each normal spring

  double *posLy; //local-y position
  double *posLz; //local-z position
  double *distFct; //distribution factor



  //orient, mass
  Vector oriX;   // local x direction
  Vector oriYp;  // local yp direction
  double mass; // mass of element

  
  
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
