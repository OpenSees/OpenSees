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
                                                                        
// $Revision: 1.5 $
// $Date: 2007/07/27 19:23:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/LehighJoint2d.h,v $
                                                                        
// Written: CY Seo
// Created August 2008
                                                                        
#ifndef LehighJoint2d_h
#define LehighJoint2d_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <FileStream.h>
#include <OPS_Stream.h>

class Node;
class Channel;
class FEM_ObjectBroker;
class Response;
class Renderer;
class UniaxialMaterial;

class LehighJoint2d: public Element
{
 public:
  // default constructor
   LehighJoint2d(); 
  
  // defined constructor
   LehighJoint2d(int tag,int Nd1, int Nd2, int Nd3, int Nd4,
		    UniaxialMaterial& theMat1, UniaxialMaterial& theMat2,
		    UniaxialMaterial& theMat3, UniaxialMaterial& theMat4,
		    UniaxialMaterial& theMat5, UniaxialMaterial& theMat6,
		    UniaxialMaterial& theMat7, UniaxialMaterial& theMat8,
			UniaxialMaterial& theMat9);
  

  // default destructor
   ~LehighJoint2d();

  
  ////////////// public methods to obtain information about dof & connectivity    
  bool	isSubdomain(void) { return false; } ;
  
  // get number of external nodes
  int getNumExternalNodes(void) const;
  
  // return connected external nodes
  const ID &getExternalNodes(void);
  Node **getNodePtrs(void);
  
  // return number of DOFs
  int getNumDOF(void);	
  
  // set domain performs check on dof and associativity with node
  void setDomain(Domain *theDomain);
  
  //////////////////////////// public methods to set the state of the element    
  
  // commit state
  int commitState(void);
  
  // revert to last commit
  int revertToLastCommit(void);        
  
  // revert to start
  int revertToStart(void);        
  
  // determine current strain and set strain in material
  int update(void);
    
  // returns converged tangent stiffness matrix
  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);   
  
  // not required for this element formulation
  const Matrix &getDamp(void);    
  const Matrix &getMass(void);    
  
  // not required for this element formulation
  void zeroLoad(void);	
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);
  
  // get converged residual
  const Vector &getResistingForce(void);
  
  // get converged residual with inertia terms
  const Vector &getResistingForceIncInertia(void);            
  
  // public methods for element output for parallel and database processing
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  
  // display element graphically
  int displaySelf(Renderer &theViewer, int displayMode, float fact);    
  
  // print out element data
  void Print(OPS_Stream &s, int flag =0);    
  
  // implemented to print into file
  const char *getClassType(void) const {return " LehighJoint2d";};

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  
  int setParameter (char **argv, int argc, Information &info);
  int updateParameter (int parameterID, Information &info);
  
  
 protected:
  
 private:
  
  // private methods
  void getBasicTrialDisp() ;
  void getAvp(void);
  
  // material info
  UniaxialMaterial **MaterialPtr;  // pointer to the 13 different materials
  
  // node info
  ID  connectedExternalNodes;   // contains the tags of the end nodes
  Node* nodePtr[4];             // pointers to four nodes
  
  
  int nodeDbTag, dofDbTag;
  int numDOF, numBasicDOF;
  
  // various other element parameters
  double elemWidth;
  double elemHeight;
  
  Vector vs;             // vector of external commited displacements
  Vector vt;          // vector of internal commited displacements   
  Matrix avp;       // matrix describing relation between the component deformations and the external and internal deformations
  Matrix apq;         // matrix of derivative of internal equilibrium 
  
  Matrix K;              // element stiffness matrix
  Vector R;              // element residual matrix  
};

#endif
