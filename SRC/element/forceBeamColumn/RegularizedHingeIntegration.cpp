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

// $Revision: 1.2 $
// $Date: 2008-12-03 23:43:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/RegularizedHingeIntegration.cpp,v $

// Theory Reference
// ----------------
// Scott, M.H. and Hamutcuoglu, O.M. "Numerically consistent regularization
// of force-based frame elements." International Journal for Numerical
// Methods in Engineering. http://dx.doi.org/10.1002/nme.2386

#include <RegularizedHingeIntegration.h>
#include <ElementalLoad.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

RegularizedHingeIntegration::RegularizedHingeIntegration(BeamIntegration &bi,
							 double lpi, double lpj,
							 double epsi, double epsj):
  BeamIntegration(BEAM_INTEGRATION_TAG_RegularizedHinge),
  beamInt(0), lpI(lpi), lpJ(lpj), epsI(epsi), epsJ(epsj), wf(0),
  parameterID(0)
{
  beamInt = bi.getCopy();
  if (beamInt == 0) {
    opserr << "RegularizedHingeIntegration::RegularizedHingeIntegration -- failed to get copy of BeamIntegration" << endln;
  }
}

RegularizedHingeIntegration::RegularizedHingeIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_RegularizedHinge),
  beamInt(0), lpI(0.0), lpJ(0.0), epsI(0.0), epsJ(0.0), wf(0), parameterID(0)
{

}

RegularizedHingeIntegration::~RegularizedHingeIntegration()
{
  if (beamInt != 0)
    delete beamInt;

  if (wf != 0)
    delete [] wf;
}

void
RegularizedHingeIntegration::getSectionLocations(int numSections, double L,
						 double *xi)
{
  beamInt->getSectionLocations(numSections-2, L, xi);
  double tmp = xi[numSections-3];

  double oneOverL = 1.0/L;

  for (int i = numSections-1; i > 3; i--)
    xi[i] = xi[i-3];
  xi[1] = epsI*oneOverL;
  xi[2] = 1.0-epsJ*oneOverL;
  xi[3] = tmp;
}

void
RegularizedHingeIntegration::getSectionWeights(int numSections, double L,
					       double *wt)
{
  beamInt->getSectionWeights(numSections-2, L, wt);

  double oneOverL = 1.0/L;

  double betaI = lpI*oneOverL;
  wt[1] = wt[0]-betaI;
  wt[0] = betaI;

  double betaJ = lpJ*oneOverL;
  wt[2] = wt[numSections-3]-betaJ;
  wt[3] = betaJ;
  
  const int nc = 4;
  int nf = numSections - nc;

  if (nf > 0) {
    
    if (wf == 0)
      wf = new double[nf];
    
    double pt[100];
    this->getSectionLocations(numSections, L, pt);

    Vector wc(wt, nc);
    Vector xc(pt, nc);
    Vector xf(&pt[nc], nf);

    Vector R(nf);
    for (int i = 0; i < nf; i++) {
      double sum = 0.0;
      for (int j = 0; j < nc; j++)
	sum += pow(xc(j),i)*wc(j);
      R(i) = 1.0/(i+1) - sum;
    }
    
    Matrix J(nf,nf);
    for (int i = 0; i < nf; i++)
      for (int j = 0; j < nf; j++)
	J(i,j) = pow(xf(j),i);
    
    Vector wfVec(wf, nf);
    
    J.Solve(R, wfVec);
  }

  for (int i = 0; i < nf; i++)
    wt[i+nc] = wf[i];
}

BeamIntegration*
RegularizedHingeIntegration::getCopy(void)
{
  return new RegularizedHingeIntegration(*beamInt, lpI, lpJ, epsI, epsJ);
}

int
RegularizedHingeIntegration::setParameter(const char **argv, int argc,
					  Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"lpI") == 0) {
    param.setValue(lpI);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"lpJ") == 0) {
    param.setValue(lpJ);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"lp") == 0) {
    param.setValue(lpI);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"zetaI") == 0) {
    param.setValue(epsI);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"zetaJ") == 0) {
    param.setValue(epsJ);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"zeta") == 0) {
    param.setValue(epsI);
    return param.addObject(6, this);
  }
  return -1;
}

int
RegularizedHingeIntegration::updateParameter(int parameterID,
					     Information &info)
{
  switch (parameterID) {
  case 1:
    lpI = info.theDouble;
    return 0;
  case 2:
    lpJ = info.theDouble;
    return 0;
  case 3:
    lpI = lpJ = info.theDouble;
    return 0;
  case 4:
    epsI = info.theDouble;
    return 0;
  case 5:
    epsJ = info.theDouble;
    return 0;
  case 6:
    epsI = epsJ = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
RegularizedHingeIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void 
RegularizedHingeIntegration::getLocationsDeriv(int numSections,
					       double L, double dLdh,
					       double *dptsdh)
{
  double oneOverL = 1.0/L;

  for (int i = 0; i < numSections; i++)
    dptsdh[i] = 0.0;

  if (parameterID == 4 || parameterID == 6) { // zetaI
    dptsdh[1] = oneOverL;
  }

  if (parameterID == 5 || parameterID == 6) { // zetaJ
    dptsdh[2] = -oneOverL;
  }

  return;

  if (dLdh != 0.0) {
    // STILL TO DO
    opserr << "getPointsDeriv -- to do" << endln;
  }

  return;
}

void
RegularizedHingeIntegration::getWeightsDeriv(int numSections,
					     double L, double dLdh,
					     double *dwtsdh)
{
  double oneOverL = 1.0/L;

  const int Nc = 4;
  int Nf = numSections - Nc;

  double dxcdh[Nc];
  double dwcdh[Nc];
  double dxfdh[100];

  for (int i = 0; i < numSections; i++) {
    dwtsdh[i] = 0.0;
    dxfdh[i] = 0.0;
  }
  for (int i = 0; i < Nc; i++) {
    dwcdh[i] = 0.0;
    dxcdh[i] = 0.0;
  }

  if (parameterID == 1 || parameterID == 3) { // lpI
    dwcdh[0] = oneOverL;
    dwcdh[1] = -oneOverL;
  }

  if (parameterID == 2 || parameterID == 3) { // lpJ
    dwcdh[2] = -oneOverL;
    dwcdh[3] = oneOverL;
  }

  if (parameterID == 4 || parameterID == 6) // epsI
    dxcdh[1] = oneOverL;

  if (parameterID == 5 || parameterID == 6) // epsJ
    dxcdh[2] = -oneOverL;

  dwtsdh[0] = dwcdh[0];
  dwtsdh[1] = dwcdh[1];
  dwtsdh[2] = dwcdh[2];
  dwtsdh[3] = dwcdh[3];

  if (Nf > 0) {

    double wt[100];
    this->getSectionWeights(numSections, L, wt);
    
    double pt[100];
    this->getSectionLocations(numSections, L, pt);

    Vector wc(wt, Nc);
    Vector xc(pt, Nc);
    Vector xf(&pt[Nc], Nf);

    Vector R(Nf);

    double sum = 0.0;
    for (int j = 0; j < Nc; j++)
      sum += dwcdh[j];
    R(0) = -sum;

    for (int i = 1; i < Nf; i++) {
      sum = 0.0;
      for (int j = 0; j < Nf; j++)
	sum += i*pow(xf(j),i-1)*dxfdh[j]*wt[Nc+j];
      for (int j = 0; j < Nc; j++)
	sum += i*pow(xc(j),i-1)*dxcdh[j]*wc[j];
      for (int j = 0; j < Nc; j++)
	sum += pow(xc(j),i)*dwcdh[j];
      R(i) = -sum;
    }

    Matrix J(Nf,Nf);
    for (int i = 0; i < Nf; i++)
      for (int j = 0; j < Nf; j++)
	J(i,j) = pow(xf(j),i);
    
    Vector dwfdh(Nf);

    J.Solve(R,dwfdh);

    for (int i = 0; i < Nf; i++)
      dwtsdh[i+Nc] = dwfdh(i);    
  }

  return;

  if (dLdh != 0.0) {
    // STILL TO DO
    opserr << "getWeightsDeriv -- to do" << endln;
  }

  return;
}

int
RegularizedHingeIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(4);

  data(0) = lpI;
  data(1) = lpJ;
  data(2) = epsI;
  data(3) = epsJ;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "RegularizedHingeIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
RegularizedHingeIntegration::recvSelf(int cTag, Channel &theChannel,
				      FEM_ObjectBroker &theBroker)
{
  static Vector data(4);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "RegularizedHingeIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  lpI = data(0);
  lpJ = data(1);
  epsI = data(2);
  epsJ = data(3);

  return 0;
}

void
RegularizedHingeIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"RegularizedHinge\", ";
		s << "\"lpI\": " << lpI << ", ";
		s << "\"lpJ\": " << lpJ << ", ";
		s << "\"epsI\": " << epsI << ", ";
		s << "\"epsJ\": " << epsJ << ", ";
		s << "\"integration\": ";
		beamInt->Print(s, flag);
		s << "}";
	}
	
	else {
		s << "RegularizedHinge" << endln;
		s << " lpI = " << lpI;
		s << " lpJ = " << lpJ << endln;
		s << " epsI = " << epsI;
		s << " epsJ = " << epsJ << endln;
		beamInt->Print(s, flag);
	}

  return;
}
