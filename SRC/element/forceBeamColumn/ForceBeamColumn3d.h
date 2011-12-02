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

// $Revision: 1.12 $
// $Date: 2010-09-13 21:26:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ForceBeamColumn3d.h,v $

/*
 * References
 *

State Determination Algorithm
---
Neuenhofer, A. and F. C. Filippou (1997). "Evaluation of Nonlinear Frame Finite
Element Models." Journal of Structural Engineering, 123(7):958-966.

Spacone, E., V. Ciampi, and F. C. Filippou (1996). "Mixed Formulation of
Nonlinear Beam Finite Element." Computers and Structures, 58(1):71-83.


Plastic Hinge Integration
---
Scott, M. H. and G. L. Fenves (2006). "Plastic Hinge Integration Methods for
Force-Based Beam-Column Elements." Journal of Structural Engineering,
132(2):244-252.


Analytical Response Sensitivity (DDM)
---
Scott, M. H., P. Franchin, G. L. Fenves, and F. C. Filippou (2004).
"Response Sensitivity for Nonlinear Beam-Column Elements."
Journal of Structural Engineering, 130(9):1281-1288.


Software Design
---
Scott, M. H., G. L. Fenves, F. T. McKenna, and F. C. Filippou (2007).
"Software Patterns for Nonlinear Beam-Column Models."
Journal of Structural Engineering, Approved for publication, February 2007.

 *
 */

#ifndef ForceBeamColumn3d_h
#define ForceBeamColumn3d_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <BeamIntegration.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>

class Response;
class ElementalLoad;

class ForceBeamColumn3d: public Element
{
 public:
  ForceBeamColumn3d();
  ForceBeamColumn3d(int tag, int nodeI, int nodeJ, 
		    int numSections, SectionForceDeformation **sec,
		    BeamIntegration &beamIntegr,
		    CrdTransf &coordTransf, double rho = 0.0, 
		    int maxNumIters = 10, double tolerance = 1.0e-12);
  
  ~ForceBeamColumn3d();
  const char *getClassType(void) const {return "ForceBeamColumn3d";};
  
  int getNumExternalNodes(void) const;
  const ID &getExternalNodes(void);
  Node **getNodePtrs(void);
  
  int getNumDOF(void);
  
  void setDomain(Domain *theDomain);
  int commitState(void);
  int revertToLastCommit(void);        
  int revertToStart(void);
  int update(void);    
  
  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getMass(void);    
  
  void zeroLoad(void);	
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);
  
  const Vector &getResistingForce(void);
  const Vector &getResistingForceIncInertia(void);            
  
  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  int displaySelf(Renderer &theViewer, int displayMode, float fact);        
  
  friend OPS_Stream &operator<<(OPS_Stream &s, ForceBeamColumn3d &E);        
  void Print(OPS_Stream &s, int flag =0);    
  
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  
 protected:
  void setSectionPointers(int numSections, SectionForceDeformation **secPtrs);
  int getInitialFlexibility(Matrix &fe);
  int getInitialDeformations(Vector &v0);
  
 private:
  void getForceInterpolatMatrix(double xi, Matrix &b, const ID &code);
  void getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code);
  void compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const;
  void initializeSectionHistoryVariables (void);
  
  // internal data
  ID     connectedExternalNodes; // tags of the end nodes

  BeamIntegration *beamIntegr;
  int numSections;
  SectionForceDeformation **sections;          // array of pointers to sections
  CrdTransf *crdTransf;        // pointer to coordinate tranformation object 
  // (performs the transformation between the global and basic system)
  double rho;                    // mass density per unit length
  int    maxIters;               // maximum number of local iterations
  double tol;	                   // tolerance for relative energy norm for local iterations
  
  int    initialFlag;            // indicates if the element has been initialized
  
  Node *theNodes[2];   // pointers to the nodes
  
  Matrix kv;                     // stiffness matrix in the basic system 
  Vector Se;                     // element resisting forces in the basic system
  
  Matrix kvcommit;               // commited stiffness matrix in the basic system
  Vector Secommit;               // commited element end forces in the basic system
  
  Matrix *fs;                    // array of section flexibility matrices
  Vector *vs;                    // array of section deformation vectors
  Vector *Ssr;                   // array of section resisting force vectors
  
  Vector *vscommit;              // array of commited section deformation vectors
  
  Matrix *sp;
  double p0[5]; // Reactions in the basic system due to element loads
  double v0[5]; // Initial deformations due to element loads

  Matrix *Ki;

  bool isTorsion;
  
  static Matrix theMatrix;
  static Vector theVector;
  static double workArea[];
  
  enum {maxNumSections = 10};
  
  // following are added for subdivision of displacement increment
  int    maxSubdivisions;       // maximum number of subdivisons of dv for local iterations
  
  static Vector *vsSubdivide;
  static Vector *SsrSubdivide;
  static Matrix *fsSubdivide;
  //static int maxNumSections;
};

#endif
