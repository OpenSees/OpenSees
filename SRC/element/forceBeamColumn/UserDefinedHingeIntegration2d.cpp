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

// $Revision: 1.3 $
// $Date: 2003-03-15 00:09:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/UserDefinedHingeIntegration2d.cpp,v $

#include <UserDefinedHingeIntegration2d.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>

UserDefinedHingeIntegration2d::UserDefinedHingeIntegration2d(int npL,
							     const Vector &ptL,
							     const Vector &wtL,
							     int npR,
							     const Vector &ptR,
							     const Vector &wtR,
							     double ee,
							     double aa,
							     double ii):
  BeamIntegration(BEAM_INTEGRATION_TAG_UserHinge2d),
  ptsL(npL), wtsL(npL), ptsR(npR), wtsR(npR),
  E(ee), A(aa), I(ii)
{
  int i;
  for (i = 0; i < npL; i++) {
    if (ptL(i) < 0.0 || ptL(i) > 1.0)
      opserr << "UserDefinedHingeIntegration2d::UserDefinedHingeIntegration2d -- point lies outside [0,1]" << endln;
    if (wtL(i) < 0.0 || wtL(i) > 1.0)
      opserr << "UserDefinedHingeIntegration2d::UserDefinedHingeIntegration2d -- weight lies outside [0,1]" << endln;
    ptsL(i) = ptL(i);
    wtsL(i) = wtL(i);
  }

  for (i = 0; i < npR; i++) {
    if (ptR(i) < 0.0 || ptR(i) > 1.0)
      opserr << "UserDefinedHingeIntegration2d::UserDefinedHingeIntegration2d -- point lies outside [0,1]" << endln;
    if (wtR(i) < 0.0 || wtR(i) > 1.0)
      opserr << "UserDefinedHingeIntegration2d::UserDefinedHingeIntegration2d -- weight lies outside [0,1]" << endln;
    ptsR(i) = ptR(i);
    wtsR(i) = wtR(i);
  }
}

UserDefinedHingeIntegration2d::UserDefinedHingeIntegration2d():
  BeamIntegration(BEAM_INTEGRATION_TAG_UserHinge2d),
  E(0.0), A(0.0), I(0.0)
{

}

UserDefinedHingeIntegration2d::~UserDefinedHingeIntegration2d()
{
  // Nothing to do
}

void
UserDefinedHingeIntegration2d::getSectionLocations(int numSections,
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
UserDefinedHingeIntegration2d::getSectionWeights(int numSections,
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
UserDefinedHingeIntegration2d::addElasticFlexibility(double L, Matrix &fElastic)
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
  double LoverEA  = Le/(E*A);
  double Lover3EI = Le/(3*E*I);
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

  fElastic(0,0) += LoverEA;
  fElastic(1,1) += ftmp(0,0);
  fElastic(1,2) += ftmp(0,1);
  fElastic(2,1) += ftmp(1,0);
  fElastic(2,2) += ftmp(1,1);

  return -1;
}

void
UserDefinedHingeIntegration2d::addElasticDeformations(ElementalLoad *theLoad,
						    double loadFactor,
						    double L, double *v0)
{
  return;
}

BeamIntegration*
UserDefinedHingeIntegration2d::getCopy(void)
{
  int npL = ptsL.Size();
  int npR = ptsR.Size();

  return new UserDefinedHingeIntegration2d(npL, ptsL, wtsL,
					 npR, ptsR, wtsR,
					 E, A, I);
}

int
UserDefinedHingeIntegration2d::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
UserDefinedHingeIntegration2d::recvSelf(int cTag, Channel &theChannel,
				      FEM_ObjectBroker &theBroker)
{
  return -1;
}

int
UserDefinedHingeIntegration2d::setParameter(const char **argv,
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
  else 
    return -1;
}

int
UserDefinedHingeIntegration2d::updateParameter(int parameterID,
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
  default:
    return -1;
  }
}

int
UserDefinedHingeIntegration2d::activateParameter(int parameterID)
{
  // For Terje to do
  return 0;
}
