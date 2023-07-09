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
// $Date: 2013-07-31 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/KikuchiBearing/KikuchiBearing.h,v $

// Written: Ken Ishii
// Created: Jan 2013
// Modified: Feb 17, 2015
//
// KikuchiBearing model
//
// Description: This file contains the class definition for KikuchiBearing.
//

#ifndef KikuchiBearing_h
#define KikuchiBearing_h

#include <Element.h>
#include <Matrix.h>

class Channel;
class UniaxialMaterial;
class Response;

class KikuchiBearing : public Element
{
 public:
  // constructor
  KikuchiBearing(int Tag, int Nd1, int Nd2,
		 int Shape, double Size, double TotalRubber, double TotalHeight,
		 int NMSS, UniaxialMaterial *MatMSS, double LimDisp,
		 int NMNS, UniaxialMaterial *MatMNS, double Lambda,
		 const Vector OriYp, const Vector OriX, double Mass,
		 bool IfPDInput, bool IfTilt,
		 double AdjCi, double AdjCj,
		 bool IfBalance, double LimFo, double LimFi, int NIter);
  
  KikuchiBearing();
  
  // destructor
  ~KikuchiBearing();
  
  // method to get class type
  const char *getClassType() const {return "KikuchiBearing";};
  
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

  UniaxialMaterial **theINodeMNSMaterials;  // nMNS^2 materials, MNS(i-Node)
  UniaxialMaterial **theJNodeMNSMaterials;  // nMNS^2 materials, MNS(j-Node)

  UniaxialMaterial **theMidMSSMaterials; // nMSS materials, MSS(mid-height)

  
  //shape
  int shape; // 1:round, 2:square
  double size; //diamiter (round shape), length of edge (square shape)
  double totalRubber; //total rubber thichness
  double totalHeight; //total height of bearing

  //MSS
  int nMSS; //number of shear springs in MSS
  double limDisp; //minimum deformation to calculate Feq and Seq (if limDisp < 0, never calculate)
  UniaxialMaterial *dmyMSSMaterial; //material to calculate Feq and Seq
  double mssFeq;//equivalent coefficient for force
  double mssSeq;//equivalent coefficient for stiffness

  double *cosTht; //arrangement of each spring (cos)
  double *sinTht; //arrangement of each spring (sin)

  double *commitDspMss; //CommitDisp


  //MNS
  int nMNS; //section is divided into (nMNS)*(nMNS) springs
  double lambda; //parameter, =(D/t)*sqrt(3G/K)
  double incA; //area of each spring

  double *posLy; //local-y position
  double *posLz; //local-z position
  double *distFct; //distribution factor

  double *commitStrnIMns; //CommitStrain
  double *commitStrnJMns;

  //stiff springs in mid-height
  double stfMidX; //stiffness, local-x direction
  double stfMidRY;//           local-ry 
  double stfMidRZ;//           local-rz 
  double stfMidRX;//           local-rx 

  double trialDspMidX; //deformation, local-x
  double trialDspMidRY;//             local-ry 
  double trialDspMidRZ;//             local-rz 
  double trialDspMidRX;//             local-rx 
  double trialFrcMidX; //force, local-x
  double trialFrcMidRY;//       local-ry 
  double trialFrcMidRZ;//       local-rz 
  double trialFrcMidRX;//      local-rx 

  double commitDspMidX; //deformation, local-x
  double commitDspMidRY;//             local-ry 
  double commitDspMidRZ;//             local-rz 
  double commitDspMidRX;//             local-rx 
  double commitFrcMidX; //force, local-x
  double commitFrcMidRY;//       local-ry 
  double commitFrcMidRZ;//       local-rz 
  double commitFrcMidRX;//       local-rx 

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

  Vector localIncrDisp;
  Vector incrDispij;
  Vector incrDispmn;
  
  Vector localForceij;

  //localDisp
  static Vector commitDij18; //i-m-n-j (18 conponents)
  static Vector trialDij18;

  static Vector commitFij;
  static Vector trialFij;

  //localStiff
  static Matrix Kij18;   //i-m-n-j 18x18 full
  static Matrix Kij18_11;//   sub  12x12
  static Matrix Kij18_12;//   sub  12x6
  static Matrix Kij18_21;//   sub  6x12
  static Matrix Kij18_22;//   sub  6x6
  static Matrix invKij18_22;//
  static Matrix Kij;     //i-j 12x12 reduced

  //localForce
  static Vector Fij;
  static Vector Fmn;
  
  //total
  static Vector stfCpnt;//stiffness (19 components)
  static Vector frcCpnt;//force (12 components)
  static Vector dspCpnt;//deformation (9 components)

  //options
  bool ifPDInput; //consider P-Delta moment (or not)
  bool ifTilt; //consider tilt of rigid link (or not)

  double adjCi; // P-Delta moment adjustment for reaction force
  double adjCj;

  bool ifBalance; //get rid of internal unbalanced force (or not)
  double limFo; //tolerance of external unbalanced force
  double limFi; //tolerance of internal unbalanced force
  int nIter; //number of iterations

  //subroutine

  //calculate Feq and Seq
  void subCalcMSSFeqSeq();

  //refer to finite displacement in the element
  void subRefFntDisp(bool ifCommit = true);

  //setTrialStrain for materials
  void subSetMaterialStrains(bool ifCommit = true);

  //calculate force components
  void subCalcFrcCpnt();

  //calculate stiffness components
  void subCalcStfCpnt();
  void subCalcStfCpntInit();
  void subCalcStfCpnt_main(bool ifInit);

  //make K18 matrix (full)
  void subMakeKij18();

  //make submatrices
  void subSubmatKij18();

  //reduct K18 matrix
  void subReductKij();

  //calculate Fij and Fmn
  void subMakeFijFmn();

  static Matrix theMatrix;
  static Vector theVector;
  static Vector theLoad;
};

#endif
