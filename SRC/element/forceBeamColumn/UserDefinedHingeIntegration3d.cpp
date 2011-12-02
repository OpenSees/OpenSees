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

// $Revision: 1.4 $
// $Date: 2003-06-10 00:36:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/UserDefinedHingeIntegration3d.cpp,v $

#include <UserDefinedHingeIntegration3d.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>

UserDefinedHingeIntegration3d::UserDefinedHingeIntegration3d(int npL,
							     const Vector &ptL,
							     const Vector &wtL,
							     int npR,
							     const Vector &ptR,
							     const Vector &wtR,
							     double ee,
							     double aa,
							     double iiz,
							     double iiy,
							     double gg,
							     double jj):
  BeamIntegration(BEAM_INTEGRATION_TAG_UserHinge3d),
  ptsL(npL), wtsL(npL), ptsR(npR), wtsR(npR),
  E(ee), A(aa), Iz(iiz), Iy(iiy), G(gg), J(jj)
{
  int i;
  for (i = 0; i < npL; i++) {
    if (ptL(i) < 0.0 || ptL(i) > 1.0)
      opserr << "UserDefinedHingeIntegration3d::UserDefinedHingeIntegration3d -- point lies outside [0,1]" << endln;
    if (wtL(i) < 0.0 || wtL(i) > 1.0)
      opserr << "UserDefinedHingeIntegration3d::UserDefinedHingeIntegration3d -- weight lies outside [0,1]" << endln;
    ptsL(i) = ptL(i);
    wtsL(i) = wtL(i);
  }

  for (i = 0; i < npR; i++) {
    if (ptR(i) < 0.0 || ptR(i) > 1.0)
      opserr << "UserDefinedHingeIntegration3d::UserDefinedHingeIntegration3d -- point lies outside [0,1]" << endln;
    if (wtR(i) < 0.0 || wtR(i) > 1.0)
      opserr << "UserDefinedHingeIntegration3d::UserDefinedHingeIntegration3d -- weight lies outside [0,1]" << endln;
    ptsR(i) = ptR(i);
    wtsR(i) = wtR(i);
  }
}

UserDefinedHingeIntegration3d::UserDefinedHingeIntegration3d():
  BeamIntegration(BEAM_INTEGRATION_TAG_UserHinge3d),
  E(0.0), A(0.0), Iz(0.0), Iy(0.0), G(0.0), J(0.0)
{

}

UserDefinedHingeIntegration3d::~UserDefinedHingeIntegration3d()
{
  // Nothing to do
}

void
UserDefinedHingeIntegration3d::getSectionLocations(int numSections,
						 double L, double *xi)
{
  int npL = ptsL.Size();
  int npR = ptsR.Size();

  int i, j;
  for (i = 0; i < npL; i++)
    xi[i] = ptsL(i);
  for (j = 0; j < npR; j++, i++)
    xi[i] = ptsR(j);
  for ( ; i < numSections; i++)
    xi[i] = 0.0;
}

void
UserDefinedHingeIntegration3d::getSectionWeights(int numSections,
					       double L, double *wt)
{
  int npL = wtsL.Size();
  int npR = wtsR.Size();

  int i, j;
  for (i = 0; i < npL; i++)
    wt[i] = wtsL(i);
  for (j = 0; j < npR; j++, i++)
    wt[i] = wtsR(j);
  for ( ; i < numSections; i++)
    wt[i] = 1.0;
}

int
UserDefinedHingeIntegration3d::addElasticFlexibility(double L, Matrix &fElastic)
{
  int npL = wtsL.Size();
  int npR = wtsR.Size();

  double betaI = 0.0;
  double betaJ = 0.0;

  int i;
  for (i = 0; i < npL; i++)
    betaI += wtsL(i);
  for (i = 0; i < npR; i++)
    betaJ += wtsR(i);

  // Length of elastic interior
  double Le = L*(1.0-betaI-betaJ);
  double Lover3EI = Le/(3*E*Iz);
  double Lover6EI = 0.5*Lover3EI;
  
  // Elastic flexibility of element interior
  static Matrix fe(2,2);
  fe(0,0) = fe(1,1) =  Lover3EI;
  fe(0,1) = fe(1,0) = -Lover6EI;
  
  // Equilibrium transformation matrix
  static Matrix B(2,2);
  B(0,0) = 1.0 - betaI;
  B(1,1) = 1.0 - betaJ;
  B(0,1) = -betaI;
  B(1,0) = -betaJ;
  
  // Transform the elastic flexibility of the element
  // interior to the basic system
  static Matrix ftmp(2,2);
  ftmp.addMatrixTripleProduct(0.0, B, fe, 1.0);

  fElastic(1,1) += ftmp(0,0);
  fElastic(1,2) += ftmp(0,1);
  fElastic(2,1) += ftmp(1,0);
  fElastic(2,2) += ftmp(1,1);

  Lover3EI = Le/(3*E*Iy);
  Lover6EI = 0.5*Lover3EI;
  fe(0,0) = fe(1,1) =  Lover3EI;
  fe(0,1) = fe(1,0) = -Lover6EI;
  ftmp.addMatrixTripleProduct(0.0, B, fe, 1.0);
  fElastic(3,3) += ftmp(0,0);
  fElastic(3,4) += ftmp(0,1);
  fElastic(4,3) += ftmp(1,0);
  fElastic(4,4) += ftmp(1,1);

  fElastic(0,0) += Le/(E*A);
  fElastic(5,5) += Le/(G*J);

  return -1;
}

void
UserDefinedHingeIntegration3d::addElasticDeformations(ElementalLoad *theLoad,
						    double loadFactor,
						    double L, double *v0)
{
  return;
}

BeamIntegration*
UserDefinedHingeIntegration3d::getCopy(void)
{
  int npL = ptsL.Size();
  int npR = ptsR.Size();

  return new UserDefinedHingeIntegration3d(npL, ptsL, wtsL,
					   npR, ptsR, wtsR,
					   E, A, Iz, Iy, G, J);
}

int
UserDefinedHingeIntegration3d::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
UserDefinedHingeIntegration3d::recvSelf(int cTag, Channel &theChannel,
				      FEM_ObjectBroker &theBroker)
{
  return -1;
}

int
UserDefinedHingeIntegration3d::setParameter(const char **argv,
					    int argc, Information &info)
{
  if (strcmp(argv[0],"E") == 0) {
    info.theType = DoubleType;
    return 1;
  }
  else if (strcmp(argv[0],"A") == 0) {
    info.theType = DoubleType;
    return 2;
  }
  else if (strcmp(argv[0],"Iz") == 0) {
    info.theType = DoubleType;
    return 3;
  }
  else if (strcmp(argv[0],"Iy") == 0) {
    info.theType = DoubleType;
    return 4;
  }
  else if (strcmp(argv[0],"G") == 0) {
    info.theType = DoubleType;
    return 5;
  }
  else if (strcmp(argv[0],"J") == 0) {
    info.theType = DoubleType;
    return 6;
  }
  else 
    return -1;
}

int
UserDefinedHingeIntegration3d::updateParameter(int parameterID,
					       Information &info)
{
  switch (parameterID) {
  case 1:
    E = info.theDouble;
    return 0;
  case 2:
    A = info.theDouble;
    return 0;
  case 3:
    Iz = info.theDouble;
    return 0;
  case 4:
    Iy = info.theDouble;
    return 0;
  case 5:
    G = info.theDouble;
    return 0;
  case 6:
    J = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
UserDefinedHingeIntegration3d::activateParameter(int parameterID)
{
  // For Terje to do
  return 0;
}

void
UserDefinedHingeIntegration3d::Print(OPS_Stream &s, int flag)
{
  s << "UserHinge3d" << endln;
  s << " E = " << E;
  s << " A = " << A;
  s << " Iz = " << Iz;
  s << " Iy = " << Iy;
  s << " G = " << G;
  s << " J = " << J << endln;
  s << " Points left hinge: " << ptsL;
  s << " Weights left hinge: " << wtsL;
  s << " Points right hinge: " << ptsR;
  s << " Weights right hinge: " << wtsR;
}
