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
// $URL$

// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ModElasticBeam3d.
// ModElasticBeam3d is a plane frame member.

#ifndef ModElasticBeam3d_h
#define ModElasticBeam3d_h

#include <Element.h>
#include <Matrix.h>
#include <Node.h>
#include <Vector.h>

class Channel;
class Information;
class CrdTransf;
class Response;
class Renderer;
class SectionForceDeformation;

class ModElasticBeam3d : public Element {
public:
  ModElasticBeam3d();
  ModElasticBeam3d(int tag, double A, double E, double G, double Jx, double Iy,
                   double Iz,
				   double k11y, double k33y, double k44y,
				   double k11z, double k33z, double k44z,
                   int Nd1, int Nd2, CrdTransf &theTransf,
                   double rho = 0.0, int cMass = 0);
  ~ModElasticBeam3d();

  const char *getClassType(void) const { return "ModElasticBeam3d"; };

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

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);
  int displaySelf(Renderer &theViewer, int displayMode, float fact,
                  const char **modes = 0, int numModes = 0);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &info);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);

private:
  double A, E, G, Jx, Iy, Iz, K11y, K33y, K44y, K11z, K33z, K44z;

  double rho;
  int cMass;

  static Matrix K;
  static Vector P;
  Vector Q;

  static Matrix kb;
  Vector q;
  double q0[5]; // Fixed end forces in basic system (no torsion)
  double p0[5]; // Reactions in basic system (no torsion)

  double wx;
  double wy;
  double wz;

  Node *theNodes[2];

  ID connectedExternalNodes;

  CrdTransf *theCoordTransf;
};

#endif
