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
//
// based on FourNodeQuad by MHS
// Written: Seweryn Kokot, Opole University of Technology, Poland
// Created: Aug 2020
//
// Description: This file contains the class definition for NineNodeQuad.

#ifndef NineNodeQuad_h
#define NineNodeQuad_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

class Node;
class NDMaterial;
class Response;

class NineNodeQuad : public Element {
public:
  NineNodeQuad(int tag, int nd1, int nd2, int nd3, int nd4, int nd5, int nd6,
               int nd7, int nd8, int nd9, NDMaterial &m, const char *type, double t,
               double pressure = 0.0, double rho = 0.0, double b1 = 0.0,
               double b2 = 0.0);
  NineNodeQuad();
  ~NineNodeQuad();

  const char *getClassType(void) const { return "NineNodeQuad"; };

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

  // public methods to obtain stiffness, mass, damping and residual information
  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getMass(void);

  void zeroLoad();
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  const Vector &getResistingForce(void);
  const Vector &getResistingForceIncInertia(void);

  // public methods for element output
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  int displaySelf(Renderer &, int mode, float fact,
                  const char **displayModes = 0, int numModes = 0);
  void Print(OPS_Stream &s, int flag = 0);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);

  int getResponse(int responseID, Information &eleInformation);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);

  // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial
  // stresses.
  friend class PyLiq1;
  friend class TzLiq1;

protected:
private:
  // private attributes - a copy for each object of the class

  static constexpr int nip = 9; // number of integration/Gauss points
  static constexpr int nnodes = 9; // number of nodes

  NDMaterial **theMaterial; // pointer to the ND material objects

  ID connectedExternalNodes; // Tags of quad nodes

  Node *theNodes[nnodes];

  static double matrixData[(nnodes*2)*(nnodes*2)]; // array data for matrix
  static Matrix K;              // Element stiffness, damping, and mass Matrix
  static Vector P;              // Element resisting force vector
  Vector Q;                     // Applied nodal loads
  double b[2];                  // Body forces

  double appliedB[2]; // Body forces applied with load pattern, C.McGann,
                      // U.Washington
  int applyLoad;      // flag for body force in load, C.McGann, U.Washington

  Vector pressureLoad; // Pressure load at nodes

  double thickness; // Element thickness
  double pressure;  // Normal surface traction (pressure) over entire element
                   // Note: positive for outward normal
  double rho;
  static double shp[3][nnodes]; // Stores shape functions and derivatives (overwritten)
  static double pts[nip][2]; // Stores quadrature points
  static double wts[nip];    // Stores quadrature weights

  // private member functions - only objects of this class can call these
  double shapeFunction(double xi, double eta);
  void setPressureLoadAtNodes(void);

  Matrix *Ki;
};

#endif
