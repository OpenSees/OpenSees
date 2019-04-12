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

// $Revision: 1.13 $
// $Date: 2009-02-05 16:28:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ForceBeamColumn2d.h,v $

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

#ifndef ForceBeamColumn2d_h
#define ForceBeamColumn2d_h

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

class ForceBeamColumn2d: public Element
{
 public:
  ForceBeamColumn2d();
  ForceBeamColumn2d(int tag, int nodeI, int nodeJ, 
		    int numSections, SectionForceDeformation **sec,
		    BeamIntegration &beamIntegr,
		    CrdTransf &coordTransf, double rho = 0.0, 
		    int maxNumIters = 10, double tolerance = 1.0e-12);
  
  ~ForceBeamColumn2d();
  
  const char *getClassType(void) const {return "ForceBeamColumn2d";};

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
  int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes=0, int numModes=0);
  
  friend OPS_Stream &operator<<(OPS_Stream &s, ForceBeamColumn2d &E);        
  void Print(OPS_Stream &s, int flag =0);    
  
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  
  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);
  const Vector &getResistingForceSensitivity(int gradNumber);
  const Matrix &getInitialStiffSensitivity(int gradNumber);
  const Matrix &getMassSensitivity(int gradNumber);
  int commitSensitivity(int gradNumber, int numGrads);
  int getResponseSensitivity(int responseID, int gradNumber,
			     Information &eleInformation);
  // AddingSensitivity:END ///////////////////////////////////////////

 protected:
  void setSectionPointers(int numSections, SectionForceDeformation **secPtrs);
  int getInitialFlexibility(Matrix &fe);
  int getInitialDeformations(Vector &v0);
  
 private:
  void getForceInterpolatMatrix(double xi, Matrix &b, const ID &code);
  void getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code);
  void compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const;
  void initializeSectionHistoryVariables (void);
  
  // Reactions of basic system due to element loads
  void computeReactions(double *p0);

  // Section forces due to element loads
  void computeSectionForces(Vector &sp, int isec);

  // internal data
  ID     connectedExternalNodes; // tags of the end nodes

  BeamIntegration *beamIntegr;
  int numSections;
  SectionForceDeformation **sections;          // array of pointers to sections
  CrdTransf *crdTransf;        // pointer to coordinate transformation object 
  // (performs the transformation between the global and basic system)
  double rho;                    // mass density per unit length
  int    maxIters;               // maximum number of local iterations
  double tol;	                   // tolerance for relative energy norm for local iterations
  
  int    initialFlag;            // indicates if the element has been initialized
  
  Node *theNodes[2];   // pointers to the nodes
  
  Matrix kv;                     // stiffness matrix in the basic system 
  Vector Se;                     // element resisting forces in the basic system
  
  Matrix kvcommit;               // committed stiffness matrix in the basic system
  Vector Secommit;               // committed element end forces in the basic system
  
  Matrix *fs;                    // array of section flexibility matrices
  Vector *vs;                    // array of section deformation vectors
  Vector *Ssr;                   // array of section resisting force vectors
  
  Vector *vscommit;              // array of committed section deformation vectors
  
  enum {maxNumEleLoads = 100};
  enum {NDM = 2};         // dimension of the problem (2d)
  enum {NND = 3};         // number of nodal dof's
  enum {NEGD = 6};         // number of element global dof's
  enum {NEBD = 3};         // number of element dof's in the basic system

  int numEleLoads; // Number of element load objects
  int sizeEleLoads;
  ElementalLoad **eleLoads;
  double *eleLoadFactors;
  Vector load;

  Matrix *Ki;
  
  static Matrix theMatrix;
  static Vector theVector;
  static double workArea[];
  
  enum {maxNumSections = 30};
  enum {maxSectionOrder = 5};

  // following are added for subdivision of displacement increment
  int    maxSubdivisions;       // maximum number of subdivisons of dv for local iterations
  
  static Vector *vsSubdivide;
  static Vector *SsrSubdivide;
  static Matrix *fsSubdivide;
  //static int maxNumSections;

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int parameterID;
  const Vector &computedqdh(int gradNumber);
  const Matrix &computedfedh(int gradNumber);
  void computeReactionSensitivity(double *dp0dh, int gradNumber);
  void computeSectionForceSensitivity(Vector &dspdh, int isec, int gradNumber);
  // AddingSensitivity:END ///////////////////////////////////////////

  Matrix tjcMass;
};

#endif
