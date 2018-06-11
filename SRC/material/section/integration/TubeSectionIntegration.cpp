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
// $Date: 2010-09-13 21:31:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/integration/TubeSectionIntegration.cpp,v $

#include <TubeSectionIntegration.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <math.h>

TubeSectionIntegration::TubeSectionIntegration(double DIAM,
					       double T,
					       int NFW,
					       int NFR):
  SectionIntegration(SECTION_INTEGRATION_TAG_Tube),
  D(DIAM), t(T), Nfwedge(NFW), Nfring(NFR), parameterID(0)
{
  
}

TubeSectionIntegration::TubeSectionIntegration():
  SectionIntegration(SECTION_INTEGRATION_TAG_Tube),
  D(0.0), t(0.0), Nfwedge(0), Nfring(0), parameterID(0)
{
  
}

TubeSectionIntegration::~TubeSectionIntegration()
{
  
}

int
TubeSectionIntegration::getNumFibers(FiberType type)
{
  return Nfwedge*Nfring;
}

int
TubeSectionIntegration::arrangeFibers(UniaxialMaterial **theMaterials,
				      UniaxialMaterial *theSteel)
{
  int numFibers = this->getNumFibers();

  for (int i = 0; i < numFibers; i++)
    theMaterials[i] = theSteel;
  
  return 0;
}

int
TubeSectionIntegration::arrangeFibers(NDMaterial **theMaterials,
				      NDMaterial *theSteel)
{
  int numFibers = this->getNumFibers();

  for (int i = 0; i < numFibers; i++)
    theMaterials[i] = theSteel;
  
  return 0;
}

void
TubeSectionIntegration::getFiberLocations(int nFibers, double *yi, double *zi)
{
  static const double pi = 3.141592653589793;
  
  double theta = pi/Nfwedge;
  double twoTheta = 2.0*theta;
  double dt = t/Nfring;

  int loc = 0;
  double rinner = 0.5*D - t;
  double Ainner = rinner*rinner*theta;
  double xinner = 2.0/3.0*rinner*sin(theta)/theta;
  for (int i = 0; i < Nfring; i++) {
    double router = rinner + (i+1)*dt;
    double Aouter = router*router*theta;
    double xouter = 2.0/3.0*router*sin(theta)/theta;
    double area = Aouter-Ainner;
    double xbar = (xouter*Aouter-xinner*Ainner)/area;
    double angle = theta;
    for (int j = 0; j < Nfwedge; j++) {
      yi[loc] = xbar*cos(angle);
      zi[loc] = xbar*sin(angle);
      angle += twoTheta;
      loc++;
    }
    Ainner = Aouter;
    xinner = xouter;    
  }
  
  return;
}

void
TubeSectionIntegration::getFiberWeights(int nFibers, double *wt)
{
  static const double pi = 3.141592653589793;
  
  double theta = pi/Nfwedge;
  double dt = t/Nfring;

  int loc = 0;
  double rinner = 0.5*D - t;
  double Ainner = rinner*rinner*theta;
  for (int i = 0; i < Nfring; i++) {
    double router = 0.5*D - t + (i+1)*dt;
    double Aouter = router*router*theta;
    double area = Aouter-Ainner;
    for (int j = 0; j < Nfwedge; j++)
      wt[loc++] = area;
    Ainner = Aouter;
  }
  
  return;
}

SectionIntegration*
TubeSectionIntegration::getCopy(void)
{
  return new TubeSectionIntegration(D, t, Nfwedge, Nfring);
}

int
TubeSectionIntegration::setParameter(const char **argv, int argc,
				     Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"d") == 0 || strcmp(argv[0],"D") == 0) {
    param.setValue(D);
    return param.addObject(1, this);    
  }
  if (strcmp(argv[0],"t") == 0) {
    param.setValue(t);
    return param.addObject(2, this);
  }

  return -1;
}

int
TubeSectionIntegration::updateParameter(int parameterID,
					Information &info)
{
  switch (parameterID) {
  case 1:
    D = info.theDouble;
    return 0;
  case 2:
    t = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
TubeSectionIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
TubeSectionIntegration::getLocationsDeriv(int nFibers, double *dyidh, double *dzidh)
{
  // Setting them both to zero for now
  for (int i = 0; i < nFibers; i++) {
    dyidh[i] = 0.0;
    dzidh[i] = 0.0;
  }
  
  return;
}

void
TubeSectionIntegration::getWeightsDeriv(int nFibers, double *dwtdh)
{
  // Setting to zero for now
  for (int i = 0; i < nFibers; i++) {
    dwtdh[i] = 0.0;
  }
  
  return;
}

void
TubeSectionIntegration::Print(OPS_Stream &s, int flag)
{
  s << "Tube" << endln;
  s << " D = "  << D;
  s << " t = " << t; 
  s << " Nfwedge = " << Nfwedge;
  s << " Nfring = " << Nfring << endln;

  return;
}

int
TubeSectionIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(4);

  data(0) = D;
  data(1) = t;
  data(2) = Nfwedge;
  data(3) = Nfring;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "TubeSectionIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
TubeSectionIntegration::recvSelf(int cTag, Channel &theChannel,
				 FEM_ObjectBroker &theBroker)
{
  static Vector data(4);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "TubeSectionIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  D   = data(0);
  t   = data(1);
  Nfwedge = (int)data(2);
  Nfring  = (int)data(3);

  return 0;
}
