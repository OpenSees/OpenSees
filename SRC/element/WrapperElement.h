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
// $Date: 2008-12-09 00:14:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/WrapperElement.h,v $

// Written: fmk                                                                         
                                                                        
#ifndef WrapperElement_h
#define WrapperElement_h

#include <Element.h>
#include <Matrix.h>

#include <elementAPI.h>

class Node;
class Channel;
class UniaxialMaterial;

class WrapperElement : public Element
{
 public:
  // constructors
  WrapperElement(const char *functName, eleObject *theEle);
  WrapperElement();

  // destructor
  virtual ~WrapperElement();
  
  // public methods for element operations
  virtual int getNumExternalNodes(void) const;
  virtual const ID &getExternalNodes(void);
  Node **getNodePtrs(void);
  
  virtual int getNumDOF(void);	
  virtual void setDomain(Domain *theDomain);
  
  virtual int commitState(void);
  virtual int revertToLastCommit(void);        
  virtual int revertToStart(void);        
  virtual int update(void);
  
  virtual const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  virtual const Matrix &getMass(void);    
  
  virtual const Vector &getResistingForce(void);
  
  // public methods for output
  virtual int sendSelf(int commitTag, Channel &theChannel);
  virtual int recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker);
  virtual int displaySelf(Renderer &theViewer, int displayMode, float fact);    
  virtual void Print(OPS_Stream &s, int flag =0);    
  
 protected:
  
 private:
  char *funcName;
  eleObject *theEle;

  Node **theNodes;       // pointer to the elements nodes
  double *u;            // to hold u^(k-1)_(n+1) as nodes don't save this
  double *R; 	          // data to hold ele resisting force
  double *K;            // data to hold ele tangent
  double *M;            // data to hold ele mass AGGIUNTA
  int eleType;          // defines which elmtnn subroutine to invoke
  Matrix *Ki;

  static Vector Rvector;    // vector to hold the applied load P
  static Matrix Kmatrix;
  static Matrix Mmatrix;	// AGGIUNTA
  static ID connectedNodes;
};

#endif




