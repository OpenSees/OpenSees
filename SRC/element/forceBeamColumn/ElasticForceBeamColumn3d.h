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

// $Revision: 1.1 $
// $Date: 2007-10-13 01:21:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ElasticForceBeamColumn3d.h,v $

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

#ifndef ElasticForceBeamColumn3d_h
#define ElasticForceBeamColumn3d_h

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

class ElasticForceBeamColumn3d: public Element
{
 public:
  ElasticForceBeamColumn3d();
  ElasticForceBeamColumn3d(int tag, int nodeI, int nodeJ, 
			   int numSections, SectionForceDeformation **sec,
			   BeamIntegration &beamIntegr,
			   CrdTransf &coordTransf, double rho = 0.0);
  
  ~ElasticForceBeamColumn3d();
  
  const char *getClassType(void) const {return "ElasticForceBeamColumn3d";};

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
  
  friend OPS_Stream &operator<<(OPS_Stream &s, ElasticForceBeamColumn3d &E);        
  void Print(OPS_Stream &s, int flag =0);    
  
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);

 protected:
  int getInitialFlexibility(Matrix &fe);
  
 private:
  void compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const;
  
  // Reactions of basic system due to element loads
  void computeReactions(double *p0);

  // Section forces due to element loads
  void computeSectionForces(Vector &sp, int isec);

  void computeBasicForces(Vector &q);

  // internal data
  ID     connectedExternalNodes; // tags of the end nodes

  enum {maxNumSections = 20};
  enum {maxSectionOrder = 10};

  BeamIntegration *beamIntegr;
  int numSections;
  SectionForceDeformation *sections[maxNumSections];          // array of pointers to sections
  CrdTransf *crdTransf;        // pointer to coordinate tranformation object 
  // (performs the transformation between the global and basic system)
  double rho;                    // mass density per unit length
  
  int    initialFlag;            // indicates if the element has been initialized
  
  Node *theNodes[2];   // pointers to the nodes
  
  Vector Se;                     // element resisting forces in the basic system
  enum {maxNumEleLoads = 100};
  enum {NDM = 3};         // dimension of the problem (2d)
  enum {NND = 6};         // number of nodal dof's
  enum {NEGD = 12};         // number of element global dof's
  enum {NEBD = 6};         // number of element dof's in the basic system

  int numEleLoads; // Number of element load objects
  ElementalLoad *eleLoads[maxNumEleLoads];

  static Matrix theMatrix;
  static Vector theVector;
  static double workArea[];
  
  int parameterID;
};

#endif
