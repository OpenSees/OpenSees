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

// $Revision: 1.8 $
// $Date: 2004-06-07 23:21:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/HingeRadauTwoBeamIntegration3d.cpp,v $

#include <HingeRadauTwoBeamIntegration3d.h>
#include <ElementalLoad.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>

HingeRadauTwoBeamIntegration3d::HingeRadauTwoBeamIntegration3d(double e,
							 double a,
							 double iz,
							 double iy,
							 double g,
							 double j,
							 double lpi,
							 double lpj):
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeRadauTwo3d),
  E(e), A(a), Iz(iz), Iy(iy), G(g), J(j), lpI(lpi), lpJ(lpj)
{
  // Nothing to do
}

HingeRadauTwoBeamIntegration3d::HingeRadauTwoBeamIntegration3d():
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeRadauTwo3d),
  E(0.0), A(0.0), Iz(0.0), Iy(0.0), G(0.0), J(0.0), lpI(0.0), lpJ(0.0)
{

}

HingeRadauTwoBeamIntegration3d::~HingeRadauTwoBeamIntegration3d()
{
  // Nothing to do
}

void
HingeRadauTwoBeamIntegration3d::getSectionLocations(int numSections, double L,
						 double *xi)
{
  double two3oneOverL = (2.0/3.0)/L;

  xi[0] = 0.0;
  xi[1] = lpI*two3oneOverL;
  xi[2] = 1.0-lpJ*two3oneOverL;
  xi[3] = 1.0;
  for (int i = 4; i < numSections; i++)
    xi[i] = 0.0;
}

void
HingeRadauTwoBeamIntegration3d::getSectionWeights(int numSections, double L,
					       double *wt)
{
  double oneOverL = 1.0/L;

  wt[0] = 0.25*lpI*oneOverL;
  wt[1] = 3.0*wt[0];
  wt[3] = 0.25*lpJ*oneOverL;
  wt[2] = 3.0*wt[3];
  for (int i = 4; i < numSections; i++)
    wt[i] = 1.0;
}

int
HingeRadauTwoBeamIntegration3d::addElasticFlexibility(double L, Matrix &fElastic)
{
  double oneOverL = 1.0/L;

  // Length of elastic interior
  double Le = L-lpI-lpJ;
  double Lover3EI = Le/(3*E*Iz);
  double Lover6EI = 0.5*Lover3EI;
  
  // Elastic flexibility of element interior
  static Matrix fe(2,2);
  fe(0,0) = fe(1,1) =  Lover3EI;
  fe(0,1) = fe(1,0) = -Lover6EI;
  
  // Equilibrium transformation matrix
  static Matrix B(2,2);
  double betaI = lpI*oneOverL;
  double betaJ = lpJ*oneOverL;
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
HingeRadauTwoBeamIntegration3d::addElasticDeformations(ElementalLoad *theLoad,
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

  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0)*loadFactor;  // Transverse
    double wz = data(1)*loadFactor;  // Transverse
    double wx = data(2)*loadFactor;  // Axial

    // Accumulate basic deformations due to uniform load on interior
    // Midpoint rule for axial
    v0[0] += wx*0.5*(L-lpI+lpJ)/(E*A)*Le;

    // Two point Gauss for bending ... will not be exact when
    // hinge lengths are not equal, but this is not a big deal!!!
    double x1 = lpI + 0.5*Le*(1.0-1.0/sqrt(3.0));
    double x2 = lpI + 0.5*Le*(1.0+1.0/sqrt(3.0));

    double Mz1 = 0.5*wy*x1*(x1-L);
    double Mz2 = 0.5*wy*x2*(x2-L);

    double My1 = -0.5*wz*x1*(x1-L);
    double My2 = -0.5*wz*x2*(x2-L);

    double Le2EIz = Le/(2*E*Iz);
    double Le2EIy = Le/(2*E*Iy);

    double b1, b2;
    b1 = x1*oneOverL;
    b2 = x2*oneOverL;
    v0[2] += Le2EIz*(b1*Mz1+b2*Mz2);
    v0[4] += Le2EIy*(b1*My1+b2*My2);

    b1 -= 1.0;
    b2 -= 1.0;
    v0[1] += Le2EIz*(b1*Mz1+b2*Mz2);
    v0[3] += Le2EIy*(b1*My1+b2*My2);
  }

  else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py = data(0)*loadFactor;
    double Pz = data(1)*loadFactor;
    double N  = data(2)*loadFactor;
    double aOverL = data(3);
    double a = aOverL*L;

    double Vy2 = Py*aOverL;
    double Vy1 = Py-Vy2;

    double Vz2 = Pz*aOverL;
    double Vz1 = Pz-Vz2;

    // Accumulate basic deformations of interior due to point load
    double M1, M2, M3;
    double b1, b2, b3;

    // Point load is on left hinge
    if (a < lpI) {
      M1 = (lpI-L)*Vy2;
      M2 = -lpJ*Vy2;

      double Le_6EI;
      Le_6EI = Le/(6*E*Iz);

      b1 = lpI*oneOverL;
      b2 = 1.0-lpJ*oneOverL;
      v0[2] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      b1 -= 1.0;
      b2 -= 1.0;
      v0[1] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      M1 = (L-lpI)*Vz2;
      M2 = lpJ*Vz2;

      Le_6EI = Le/(6*E*Iy);

      b1 = lpI*oneOverL;
      b2 = 1.0-lpJ*oneOverL;
      v0[4] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      b1 -= 1.0;
      b2 -= 1.0;
      v0[3] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      // Nothing to do for axial
      //v0[0] += 0.0;
    }
    // Point load is on right hinge
    else if (a > L-lpJ) {
      M1 = -lpI*Vy1;
      M2 = (lpJ-L)*Vy1;

      double Le_6EI;
      Le_6EI = Le/(6*E*Iz);

      b1 = lpI*oneOverL;
      b2 = 1.0-lpJ*oneOverL;
      v0[2] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      b1 -= 1.0;
      b2 -= 1.0;
      v0[1] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      M1 = lpI*Vz1;
      M2 = (L-lpJ)*Vz1;

      Le_6EI = Le/(6*E*Iy);

      b1 = lpI*oneOverL;
      b2 = 1.0-lpJ*oneOverL;
      v0[4] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      b1 -= 1.0;
      b2 -= 1.0;
      v0[3] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      v0[0] += N*Le/(E*A);      
    }
    // Point load is on elastic interior
    else {
      M1 = -lpI*Vy1;
      M2 = -lpJ*Vy2;
      M3 = -a*Vy1;

      double L1_6EI;
      double L2_6EI;

      L1_6EI = (a-lpI)/(6*E*Iz);
      L2_6EI = (Le-a+lpI)/(6*E*Iz);

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

      M1 = lpI*Vz1;
      M2 = lpJ*Vz2;
      M3 = a*Vz1;

      L1_6EI = (a-lpI)/(6*E*Iy);
      L2_6EI = (Le-a+lpI)/(6*E*Iy);

      b1 = lpI*oneOverL;
      b2 = 1.0-lpJ*oneOverL;
      b3 = a*oneOverL;
      v0[4] += L1_6EI*(M1*(2*b1+b3)+M3*(b1+2*b3));
      v0[4] += L2_6EI*(M2*(2*b2+b3)+M3*(b2+2*b3));

      b1 -= 1.0;
      b2 -= 1.0;
      b3 -= 1.0;
      v0[3] += L1_6EI*(M1*(2*b1+b3)+M3*(b1+2*b3));
      v0[3] += L2_6EI*(M2*(2*b2+b3)+M3*(b2+2*b3));

      v0[0] += N*(a-lpI)/(E*A);
    }
  }

  return;
}

double
HingeRadauTwoBeamIntegration3d::getTangentDriftI(double L, double LI,
						 double q2, double q3, bool yAxis)
{
  double oneOverL = 1.0/L;

  double betaI = lpI*oneOverL;

  double qq2 = (1-betaI)*q2 - betaI*q3;

  if (LI < lpI)
    return 0.0;
  else {
    double EI = yAxis ? E*Iy : E*Iz;
    return (LI-lpI)/3*(LI-lpI)*qq2/(EI);
  }
}

double
HingeRadauTwoBeamIntegration3d::getTangentDriftJ(double L, double LI,
						 double q2, double q3, bool yAxis)
{
  double oneOverL = 1.0/L;

  double betaJ = lpJ*oneOverL;

  double qq3 = (1-betaJ)*q3 - betaJ*q2;

  if (LI > L-lpJ)
    return 0.0;
  else {
    double EI = yAxis ? E*Iy : E*Iz;
    return (L-LI-lpJ)/3*(L-LI-lpJ)*qq3/EI;
  }
}

BeamIntegration*
HingeRadauTwoBeamIntegration3d::getCopy(void)
{
  return new HingeRadauTwoBeamIntegration3d(E, A, Iz, Iy, G, J, lpI, lpJ);
}

int
HingeRadauTwoBeamIntegration3d::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(8);

  data(0) = E;
  data(1) = A;
  data(2) = Iz;
  data(3) = Iy;
  data(4) = G;
  data(5) = J;
  data(6) = lpI;
  data(7) = lpJ;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "HingeRadauTwoBeamIntegration3d::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
HingeRadauTwoBeamIntegration3d::recvSelf(int cTag, Channel &theChannel,
				      FEM_ObjectBroker &theBroker)
{
  static Vector data(8);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "HingeRadauTwoBeamIntegration3d::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  E   = data(0);
  A   = data(1);
  Iz  = data(2);
  Iy  = data(3);
  G   = data(4);
  J   = data(5);
  lpI = data(6);
  lpJ = data(7);

  return 0;
}

int
HingeRadauTwoBeamIntegration3d::setParameter(const char **argv,
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
  else if (strcmp(argv[0],"lpI") == 0) {
    info.theType = DoubleType;
    return 7;
  }
  else if (strcmp(argv[0],"lpJ") == 0) {
    info.theType = DoubleType;
    return 8;
  }
  else 
    return -1;
}

int
HingeRadauTwoBeamIntegration3d::updateParameter(int parameterID,
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
  case 7:
    lpI = info.theDouble;
    return 0;
  case 8:
    lpJ = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
HingeRadauTwoBeamIntegration3d::activateParameter(int parameterID)
{
  // For Terje to do
  return 0;
}

void
HingeRadauTwoBeamIntegration3d::Print(OPS_Stream &s, int flag)
{
  s << "HingeRadauTwo3d" << endln;
  s << " E = " << E;
  s << " A = " << A;
  s << " Iz = " << Iz;
  s << " Iy = " << Iy;
  s << " G = " << G;
  s << " J = " << J << endln;
  s << " lpI = " << lpI;
  s << " lpJ = " << lpJ << endln;

  return;
}
