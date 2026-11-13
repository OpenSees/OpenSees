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

// Written: Nobuyuki Takimoto
//
// Three node flat shell element for thin and moderately thick shells.
// Ref: Plate bending part - DKMT (Discrete Kirchhoff-Mindlin Triangle),
//      I. Katili, "A new discrete Kirchhoff-Mindlin element based on
//      Mindlin-Reissner plate theory and assumed shear strain fields -
//      Part I: An extended DKT element for thick-plate bending analysis",
//      IJNME 36 (1993) 1859-1883. Reduces to DKT when the plate is thin.
//      Membrane part - CST (constant strain triangle). The drilling DOF is
//      tied to the in-plane rigid rotation of the CST by a small penalty.

#ifndef ShellDKMT_h
#define ShellDKMT_h

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>

class ShellDKMT : public Element {

 public:

  ShellDKMT();
  ShellDKMT(int tag, int node1, int node2, int node3,
            SectionForceDeformation &theSection, double drillRatio = 1.0e-4);
  virtual ~ShellDKMT();

  const char *getClassType(void) const {return "ShellDKMT";}

  void setDomain(Domain *theDomain);

  int getNumExternalNodes() const;
  const ID &getExternalNodes();
  Node **getNodePtrs();
  int getNumDOF();

  int commitState();
  int revertToLastCommit();
  int revertToStart();

  void Print(OPS_Stream &s, int flag);

  const Matrix &getTangentStiff();
  const Matrix &getInitialStiff();
  const Matrix &getMass();

  void zeroLoad(void);
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  const Vector &getResistingForce();
  const Vector &getResistingForceIncInertia();

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  Response *setResponse(const char **argv, int argc, OPS_Stream &output);
  int getResponse(int responseID, Information &eleInfo);

  int setParameter(const char **argv, int argc, Parameter &param);

  int displaySelf(Renderer &, int mode, float fact, const char **displayModes = 0, int numModes = 0);

 private:

  static const int numGauss = 3;

  static Matrix stiff;
  static Vector resid;
  static Matrix mass;

  ID connectedExternalNodes;
  Node *nodePointers[3];
  SectionForceDeformation *sections[numGauss];

  double drillRatio;   // drilling stiffness = drillRatio * D(2,2) of the section

  // local frame (rows of the rotation matrix) and local nodal coordinates
  double g[3][3];
  double xl[3][2];
  double area;

  // strain-displacement matrices in local coordinates, 8 x 18 at each
  // Gauss point. B gives the generalized strains as expected by the shell
  // sections, B1 is B with the bending rows inverted (the sign convention
  // of the shell sections, as in ASDShellT3)
  Matrix B[numGauss];
  Matrix B1[numGauss];
  double Bdrill[3][18];   // theta_z_i - in-plane rigid rotation of the CST
  double kdrill;

  Vector *load;
  double appliedB[3];
  int applyLoad;

  void computeBasis();
  void formB();
  void formLocal(Matrix &K, Vector &R, int tangFlag);
  void formInertiaTerms(int tangFlag);
  void localToGlobal(const Matrix &Kl, const Vector &Rl, Matrix &K, Vector &R, bool doK);
  void localDisp(Vector &ul);
};

#endif
