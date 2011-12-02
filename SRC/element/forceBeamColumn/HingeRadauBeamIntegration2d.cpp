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

// $Revision: 1.10 $
// $Date: 2004-06-07 23:21:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/HingeRadauBeamIntegration2d.cpp,v $
#include <HingeRadauBeamIntegration2d.h>
#include <ElementalLoad.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>

HingeRadauBeamIntegration2d::HingeRadauBeamIntegration2d(double e,
							 double a,
							 double i,
							 double lpi,
							 double lpj):
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeRadau2d),
  E(e), A(a), I(i), lpI(lpi), lpJ(lpj)
{
  // Nothing to do
}

HingeRadauBeamIntegration2d::HingeRadauBeamIntegration2d():
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeRadau2d),
  E(0.0), A(0.0), I(0.0), lpI(0.0), lpJ(0.0)
{

}

HingeRadauBeamIntegration2d::~HingeRadauBeamIntegration2d()
{
  // Nothing to do
}

void
HingeRadauBeamIntegration2d::getSectionLocations(int numSections, double L,
						 double *xi)
{
  xi[0] = 0.0;
  xi[1] = 1.0;
  for (int i = 2; i < numSections; i++)
    xi[i] = 0.0;
}

void
HingeRadauBeamIntegration2d::getSectionWeights(int numSections, double L,
					       double *wt)
{
  double oneOverL = 1.0/L;

  wt[0] = lpI*oneOverL;
  wt[1] = lpJ*oneOverL;
  for (int i = 2; i < numSections; i++)
    wt[i] = 1.0;
}

int
HingeRadauBeamIntegration2d::addElasticFlexibility(double L, Matrix &fElastic)
{
  double oneOverL = 1.0/L;

  double Le = L-lpI-lpJ;
  fElastic(0,0) += Le/(E*A);

  double x[4];
  double w[4];
  
  static const double eight3 = 8.0/3.0;
  static const double oneOverRoot3 = 1.0/sqrt(3.0);

  double oneOverEI = 1.0/(E*I);
  
  x[0] = eight3*lpI;
  w[0] = 3.0*lpI;

  x[1] = L-eight3*lpJ;
  w[1] = 3.0*lpJ;

  Le = L-4.0*(lpI+lpJ);

  x[2] = 4.0*lpI + 0.5*Le*(1.0-oneOverRoot3);
  w[2] = 0.5*Le;

  x[3] = 4.0*lpI + 0.5*Le*(1.0+oneOverRoot3);
  w[3] = w[2];

  double tmp = 0.0;
  double xL, xL1, wt;
  for (int i = 0; i < 4; i++) {
    xL  = x[i]*oneOverL;
    xL1 = xL-1.0;
    wt = w[i]*oneOverEI;
    fElastic(1,1) += xL1*xL1*wt;
    fElastic(2,2) += xL*xL*wt;
    tmp           += xL*xL1*wt;
  }
  fElastic(1,2) += tmp;
  fElastic(2,1) += tmp;

  return -1;
}

void
HingeRadauBeamIntegration2d::addElasticDeformations(ElementalLoad *theLoad,
						    double loadFactor,
						    double L, double *v0)
{
  // Length of elastic interior
  double Le = L-lpI-lpJ;
  if (Le == 0.0)
    return;

  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  double oneOverL = 1.0/L;

  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wa = data(1)*loadFactor;  // Axial
    double wy = data(0)*loadFactor;  // Transverse

    // Accumulate basic deformations due to uniform load on interior
    // Midpoint rule for axial
    v0[0] += wa*0.5*(L-lpI+lpJ)/(E*A)*Le;

    // Two point Gauss for bending ... will not be exact when
    // hinge lengths are not equal, but this is not a big deal!!!
    double x1 = lpI + 0.5*Le*(1.0-1.0/sqrt(3.0));
    double x2 = lpI + 0.5*Le*(1.0+1.0/sqrt(3.0));

    double M1 = 0.5*wy*x1*(x1-L);
    double M2 = 0.5*wy*x2*(x2-L);

    double b1, b2;
    double Le2EI = Le/(2*E*I);
    b1 = x1*oneOverL;
    b2 = x2*oneOverL;
    v0[2] += Le2EI*(b1*M1+b2*M2);

    b1 -= 1.0;
    b2 -= 1.0;
    v0[1] += Le2EI*(b1*M1+b2*M2);
  }

  else if (type == LOAD_TAG_Beam2dPointLoad) {
    double P = data(0)*loadFactor;
    double N = data(1)*loadFactor;
    double aOverL = data(2);
    double a = aOverL*L;

    double V1 = P*(1.0-aOverL);
    double V2 = P*aOverL;

    // Accumulate basic deformations of interior due to point load
    double M1, M2, M3;
    double b1, b2, b3;

    // Point load is on left hinge
    if (a < lpI) {
      M1 = (lpI-L)*V2;
      M2 = -lpJ*V2;

      double Le_6EI = Le/(6*E*I);

      b1 = lpI*oneOverL;
      b2 = 1.0-lpJ*oneOverL;
      v0[2] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      b1 -= 1.0;
      b2 -= 1.0;
      v0[1] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      // Nothing to do for axial
      //v0[0] += 0.0;
    }
    // Point load is on right hinge
    else if (a > L-lpJ) {
      M1 = -lpI*V1;
      M2 = (lpJ-L)*V1;

      double Le_6EI = Le/(6*E*I);

      b1 = lpI*oneOverL;
      b2 = 1.0-lpJ*oneOverL;
      v0[2] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      b1 -= 1.0;
      b2 -= 1.0;
      v0[1] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      v0[0] += N*Le/(E*A);      
    }
    // Point load is on elastic interior
    else {
      M1 = -lpI*V1;
      M2 = -lpJ*V2;
      M3 = -a*V1;

      double L1_6EI = (a-lpI)/(6*E*I);
      double L2_6EI = (Le-a+lpI)/(6*E*I);

      b1 = lpI*oneOverL;
      b2 = 1.0-lpJ*oneOverL;
      b3 = a*oneOverL;
      v0[2] += L1_6EI*(M1*(2*b1+b3)+M3*(b1+2*b3));
      v0[2] += L2_6EI*(M2*(2*b2+b3)+M3*(b2+2*b3));

      b1 -= 1.0;
      b2 -= 1.0;
      b3 -= 1.0;
      v0[1] += L1_6EI*(M1*(2*b1+b3)+M3*(b1+2*b3));
      v0[1] += L2_6EI*(M2*(2*b2+b3)+M3*(b2+2*b3));

      v0[0] += N*(a-lpI)/(E*A);
    }
  }

  return;
}

double
HingeRadauBeamIntegration2d::getTangentDriftI(double L, double LI,
					      double q2, double q3)
{
  double oneOverL = 1.0/L;

  double betaI = 4*lpI*oneOverL;

  double qq2 = (1-betaI)*q2 - betaI*q3;

  betaI = 8.0/3*lpI*oneOverL;

  double qqq2 = (1-betaI)*q2 - betaI*q3;

  if (LI < lpI)
    return 0.0;
  else
    return (3*lpI)*(LI-8.0/3*lpI)*qqq2/(E*I) +
      (LI-4*lpI)/3*(LI-4*lpI)*qq2/(E*I);
}

double
HingeRadauBeamIntegration2d::getTangentDriftJ(double L, double LI,
					      double q2, double q3)
{
  double oneOverL = 1.0/L;

  double betaJ = 4*lpJ*oneOverL;

  double qq3 = (1-betaJ)*q3 - betaJ*q2;

  betaJ = 8.0/3*lpJ*oneOverL;

  double qqq3 = (1-betaJ)*q3 - betaJ*q2;

  if (LI > L-lpJ)
    return 0.0;
  else
    return (3*lpJ)*(L-LI-8.0/3*lpJ)*qqq3/(E*I) +
      (L-LI-4*lpJ)/3*(L-LI-4*lpJ)*qq3/(E*I);
}

BeamIntegration*
HingeRadauBeamIntegration2d::getCopy(void)
{
  return new HingeRadauBeamIntegration2d(E, A, I, lpI, lpJ);
}

int
HingeRadauBeamIntegration2d::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(5);

  data(0) = E;
  data(1) = A;
  data(2) = I;
  data(3) = lpI;
  data(4) = lpJ;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "HingeRadauBeamIntegration2d::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
HingeRadauBeamIntegration2d::recvSelf(int cTag, Channel &theChannel,
				      FEM_ObjectBroker &theBroker)
{
  static Vector data(5);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "HingeRadauBeamIntegration2d::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  E   = data(0);
  A   = data(1);
  I   = data(2);
  lpI = data(3);
  lpJ = data(4);

  return 0;
}

int
HingeRadauBeamIntegration2d::setParameter(const char **argv,
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
  else if (strcmp(argv[0],"I") == 0 || strcmp(argv[0],"Iz") == 0) {
    info.theType = DoubleType;
    return 3;
  }
  else if (strcmp(argv[0],"lpI") == 0) {
    info.theType = DoubleType;
    return 4;
  }
  else if (strcmp(argv[0],"lpJ") == 0) {
    info.theType = DoubleType;
    return 5;
  }
  else 
    return -1;
}

int
HingeRadauBeamIntegration2d::updateParameter(int parameterID,
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
    I = info.theDouble;
    return 0;
  case 4:
    lpI = info.theDouble;
    return 0;
  case 5:
    lpJ = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
HingeRadauBeamIntegration2d::activateParameter(int parameterID)
{
  // For Terje to do
  return 0;
}

void
HingeRadauBeamIntegration2d::Print(OPS_Stream &s, int flag)
{
  s << "HingeRadau2d" << endln;
  s << " E = " << E;
  s << " A = " << A;
  s << " I = " << I << endln;
  s << " lpI = " << lpI;
  s << " lpJ = " << lpJ << endln;

  return;
}
