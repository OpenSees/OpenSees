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
                                                                        
#ifndef MultiFP2d_h
#define MultiFP2d_h

// Written: fmk
// Conversion from matlab to c++: fmk
//
// Description: This file contains the interface for the MultiFP2d class.
//
// What: "@(#) MultiFP2d.h, revA"

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <UniaxialMaterial.h>

class MultiFP2d : public Element
{
  public:
    // constructors
  MultiFP2d(int tag, 
	    int Nd1, int Nd2, 
	    UniaxialMaterial *theFrictionModel,
	    UniaxialMaterial *theVerticalModel,
	    double w0, int axialCase);

  MultiFP2d(int tag, 
	    int Nd1, int Nd2,
	    int type,
	    const Vector &R,
	    const Vector &h, 
	    const Vector &D,
	    const Vector &d,
	    const Vector &mu,
	    double Kvert,
	    double w0, int axialCase);
  
  MultiFP2d();    
  
  // destructor
  ~MultiFP2d();
  
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
  int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);  

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  
 protected:
  
 private:
  UniaxialMaterial *theFrictionModel;
  UniaxialMaterial *theVerticalModel;
  
  double W0;   // original weight
  double cW;   // current vertical force (obtained form vertical spring)
  
  ID  externalNodes;  // contains the id's of end nodes
  Node *theNodes[2];  // node pointers
  
  int numDOF;
  Matrix *theMatrix;
  Vector *theVector;  

  int type;
  int axialCase;
  Matrix data;
};
#endif

