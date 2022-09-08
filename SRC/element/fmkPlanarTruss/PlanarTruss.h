#ifndef PlanarTruss_h
#define PlanarTruss_h

/********************************************************************************
(C) Copyright 2001-2022, The Regents of the University of California    
All Rights Reserved.                                               

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or other
materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.

******************************************************************************* */


// Written: fmk 
// Created: 08/01

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class UniaxialMaterial;

#define ELE_TAG_PlanarTruss 100001
#include <UniaxialMaterial.h>

class PlanarTruss : public Element
{
  public:
  PlanarTruss(int tag, int iNode, int jNode, UniaxialMaterial *theMat, double a);
    PlanarTruss();    
    ~PlanarTruss();

    const char *getClassType(void) const {return "PlanarTruss";};

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
    
    // public methods to obtain stiffness, mass, damping and 
    // residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);    

    const Vector &getResistingForce(void);

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &);
    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:
  ID  connectedExternalNodes;   // contains the tags of the end nodes
  Matrix theMatrix;             // matrix to return stiff, damp & mass
  Vector theVector;             // vector to return the residual
  UniaxialMaterial *theMaterial;
  
  double A, L;
  Matrix trans;
  Node *theNodes[2];
};

#endif




