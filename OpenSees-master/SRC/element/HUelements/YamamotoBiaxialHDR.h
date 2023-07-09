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
// $Date: 2013-04-25 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/yamamotoBiaxialHDR/YamamotoBiaxialHDR.h,v $

// Written:  Masaru Kikuchi
// Created:  August 2014
//
// Description: This file contains the class definition for YamamotoBiaxialHDR.
//

#ifndef YamamotoBiaxialHDR_h
#define YamamotoBiaxialHDR_h

#include <Element.h>
#include <Matrix.h>

class Channel;
class Response;

class YamamotoBiaxialHDR : public Element
{
 public:
  // constructor
  YamamotoBiaxialHDR(int Tag, int Nd1, int Nd2, int Tp, double DDo, double DDi, double Hr,
	  double Cr, double Cs, const Vector OriYp, const Vector OriX = 0,
	  double Mass = 0.0);

  YamamotoBiaxialHDR();
  
  // destructor
  ~YamamotoBiaxialHDR();
  
  // method to get class type
  const char *getClassType() const {return "YamamotoBiaxialHDR";};
  
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
  //-------------------------------------------------------------------------------
  int setTrialStrain(const Vector &strain);
  const double &getStrain(int direction);
  const double &getStress(int direction);
  const double &getTangent(int direction);
  const double &getInitialTangent(int direction);
  //-------------------------------------------------------------------------------

  // private methods
  void setUp();
  
  // private attributes - a copy for each object of the class
  ID connectedExternalNodes;        // contains the tags of the end nodes
  Node *theNodes[2];                // array of nodes
  //  UniaxialMaterial **theMaterials;  // material‚ðnSpringŒÂ

  // parameters
  //  int nSpring=1;
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

  // input data
  int tp; // =1; Bridgestone X0.6R
  double ddo; // outer diameter [m]
  double ddi; // inner diameter [m]
  double hr; // total thickness of rubber layer [m]

  // model parameters
  double ar;    // cross-section area [m^2]
  double ip;    // polar moment of inertia of area [m^4]
  double alpha; // Yamamoto model parameter
  double nn; // Yamamoto model parameter
  double cr; // coefficient of tau-r
  double cs; // coefficient of tau-s

  double initialStiff[2];

  // trial values
  double trialStiff[2];
  double trialDeform[2];
  double trialForce[2];
  double trialQ[2];
  double trialP[2];
  double trialFr[2];
  double trialFs[2];


  // commit values
  double commitStiff[2];
  double commitDeform[2];
  double commitForce[2];
  double commitQ[2];
  double commitP[2];
  double commitFr[2];
  double commitFs[2];

};

#endif
