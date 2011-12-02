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
// $Date: 2008-12-03 23:42:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/LowOrderBeamIntegration.cpp,v $

#include <LowOrderBeamIntegration.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

LowOrderBeamIntegration::LowOrderBeamIntegration(int nIP,
						 const Vector &pt,
						 int nc,
						 const Vector &wc):
  BeamIntegration(BEAM_INTEGRATION_TAG_LowOrder),
  pts(nIP), wts(nIP), Nc(nc), parameterID(0)
{
  for (int i = 0; i < nIP; i++) {
    if (pt(i) < 0.0 || pt(i) > 1.0)
      opserr << "LowOrderBeamIntegration::LowOrderBeamIntegration -- point lies outside [0,1]" << endln;

    pts(i) = pt(i);
  }

  int nf = nIP-nc;

  if (nf > 0) {
    Vector R(nf);
    for (int i = 0; i < nf; i++) {
      double sum = 0.0;
      for (int j = 0; j < nc; j++)
	sum += pow(pts(j),i)*wc(j);
      R(i) = 1.0/(i+1) - sum;
    }
    
    Matrix J(nf,nf);
    for (int i = 0; i < nf; i++)
      for (int j = 0; j < nf; j++)
	J(i,j) = pow(pts(nc+j),i);
    
    Vector wf(nf);
    
    J.Solve(R, wf);
    
    for (int i = 0; i < nf; i++)
      wts(nc+i) = wf(i);
    
    for (int i = 0; i < nc; i++)
      wts(i) = wc(i);
  }
  else
    wts = wc;
}

LowOrderBeamIntegration::LowOrderBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_LowOrder)
{
 
}

LowOrderBeamIntegration::~LowOrderBeamIntegration()
{
  // Nothing to do
}

void
LowOrderBeamIntegration::getSectionLocations(int numSections,
					     double L, double *xi)
{
  int nIP = pts.Size();

  int i;
  for (i = 0; i < nIP; i++)
    xi[i] = pts(i);
  for ( ; i < numSections; i++)
    xi[i] = 0.0;
}

void
LowOrderBeamIntegration::getSectionWeights(int numSections,
					   double L, double *wt)
{
  int nIP = wts.Size();

  int Nf = nIP-Nc;

  if (Nf > 0) {
    Vector R(Nf);
    for (int i = 0; i < Nf; i++) {
      double sum = 0.0;
      for (int j = 0; j < Nc; j++)
	sum += pow(pts(j),i)*wts(j);
      R(i) = 1.0/(i+1) - sum;
    }
    
    Matrix J(Nf,Nf);
    for (int i = 0; i < Nf; i++)
      for (int j = 0; j < Nf; j++)
	J(i,j) = pow(pts(Nc+j),i);
    
    Vector wf(Nf);
    
    J.Solve(R, wf);
    
    for (int i = 0; i < Nf; i++)
      wts(Nc+i) = wf(i);
  }

  int i;
  for (i = 0; i < nIP; i++)
    wt[i] = wts(i);
  for ( ; i < numSections; i++)
    wt[i] = 1.0;
}

BeamIntegration*
LowOrderBeamIntegration::getCopy(void)
{
  int nIP = pts.Size();

  return new LowOrderBeamIntegration(nIP, pts, Nc, wts);
}

int
LowOrderBeamIntegration::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
LowOrderBeamIntegration::recvSelf(int cTag, Channel &theChannel,
				  FEM_ObjectBroker &theBroker)
{
  return -1;
}

int
LowOrderBeamIntegration::setParameter(const char **argv, int argc,
				      Parameter &param)
{
  if (argc < 2)
    return -1;

  int point = atoi(argv[1]);
  if (point < 1)
    return -1;

  int N = pts.Size();
  int Nf = N-Nc;

  if (strcmp(argv[0],"xf") == 0 && point <= Nf)
    return param.addObject(point, this);

  else if (strcmp(argv[0],"xc") == 0 && point <= Nc)
    return param.addObject(10+point, this);

  else if (strcmp(argv[0],"wc") == 0 && point <= Nc)
    return param.addObject(20+point, this);

  else
    return -1;
}

int
LowOrderBeamIntegration::updateParameter(int parameterID,
					 Information &info)
{
  if (parameterID < 10) { // xf
    pts(Nc+(parameterID-1)) = info.theDouble;
    return 0;
  }
  else if (parameterID < 20) { // xc
    pts(parameterID-10-1) = info.theDouble;
    return 0;
  }
  else if (parameterID < 30) { // wc
    pts(parameterID-20-1) = info.theDouble;
    return 0;
  }
  else
    return -1;
}

int
LowOrderBeamIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
LowOrderBeamIntegration::Print(OPS_Stream &s, int flag)
{
  s << "LowOrder" << endln;
  s << " Points: " << pts;
  s << " Weights: " << wts;
  double sum = 0.0;
  int N = wts.Size();
  for (int i = 0; i < N; i++)
    sum += fabs(wts(i));
  s << " Condition Number: " << sum << endln;
}

void 
LowOrderBeamIntegration::getLocationsDeriv(int numSections, double L,
					   double dLdh, double *dptsdh)
{
  for (int i = 0; i < numSections; i++)
    dptsdh[i] = 0.0;

  if (parameterID == 0)
    return;

  double oneOverL = 1.0/L;

  if (parameterID < 10) // xf
    dptsdh[Nc+(parameterID-1)] = oneOverL;
  else if (parameterID < 20) // xc
    dptsdh[parameterID-10-1] = oneOverL;

  return;
}

void
LowOrderBeamIntegration::getWeightsDeriv(int numSections, double L,
					 double dLdh, double *dwtsdh)
{
  for (int i = 0; i < numSections; i++)
    dwtsdh[i] = 0.0;

  if (parameterID == 0)
    return;

  double dxcdh[10];
  double dxfdh[10];
  for (int i = 0; i < 10; i++) {
    dxcdh[i] = 0.0;
    dxfdh[i] = 0.0;
  }
  
  double oneOverL = 1.0/L;

  if (parameterID < 10) // xf
    dxfdh[parameterID-1] = oneOverL;
  else if (parameterID < 20) // xc
    dxcdh[parameterID-10-1] = oneOverL;
  else if (parameterID < 30) // wc
    dwtsdh[parameterID-20-1] = oneOverL;

  int N = pts.Size();
  int Nf = N-Nc;

  if (Nf > 0) {

    Vector R(Nf);

    double sum = 0.0;
    for (int j = 0; j < Nc; j++)
      sum += dwtsdh[j];
    R(0) = -sum;

    for (int i = 1; i < Nf; i++) {
      sum = 0.0;
      for (int j = 0; j < Nf; j++)
	sum += i*pow(pts(Nc+j),i-1)*dxfdh[j]*wts(Nc+j);
      for (int j = 0; j < Nc; j++)
	sum += i*pow(pts(j),i-1)*dxcdh[j]*wts(j);
      for (int j = 0; j < Nc; j++)
	sum += pow(pts(j),i)*dwtsdh[j];
      R(i) = -sum;
    }

    Matrix J(Nf,Nf);
    for (int i = 0; i < Nf; i++)
      for (int j = 0; j < Nf; j++)
	J(i,j) = pow(pts(Nc+j),i);
    
    Vector dwfdh(Nf);

    J.Solve(R,dwfdh);

    for (int i = 0; i < Nf; i++)
      dwtsdh[Nc+i] += dwfdh(i);    
  }

  return;
}
