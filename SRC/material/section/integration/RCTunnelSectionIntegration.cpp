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
// $Source$

#include <RCTunnelSectionIntegration.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <math.h>

RCTunnelSectionIntegration::RCTunnelSectionIntegration(double D, double H,
							   double ASin, double ASout,
							   double COVERin, double COVERout,
							   int NRINGS, int NWEDGES,
							   int NBARin, int NBARout):
  SectionIntegration(SECTION_INTEGRATION_TAG_RCTUNNEL),
  d(D), h(H), Asinner(ASin), Asouter(ASout), coverinner(COVERin), coverouter(COVERout),
  Nrings(NRINGS), Nwedges(NWEDGES), Nbarsinner(NBARin), Nbarsouter(NBARout),
  parameterID(0)
{
  /*
  if (Nrings < 1)
    Nrings = 1;

  if (Nwedges < 2)
    Nwedges = 2;

  if (Nbarsinner < 1)
    Nbarsinner = 1;

  if (Nbarsouter < 1)
    Nbarsouter = 1;
  */
}

RCTunnelSectionIntegration::RCTunnelSectionIntegration():
  SectionIntegration(SECTION_INTEGRATION_TAG_RCTUNNEL),
  d(0.0), h(0.0), Asinner(0.0), Asouter(0.0), coverinner(0.0), coverouter(0.0),
  Nrings(1), Nwedges(2), Nbarsinner(1), Nbarsouter(1),
  parameterID(0)
{
  
}

RCTunnelSectionIntegration::~RCTunnelSectionIntegration()
{
  
}

int
RCTunnelSectionIntegration::getNumFibers(FiberType type)
{
  if (type == steel)
    return Nbarsinner + Nbarsouter;
  if (type == concrete)
    return Nrings*Nwedges;
  if (type == all)
    return Nrings*Nwedges + Nbarsinner + Nbarsouter;

  return 0;
}

int
RCTunnelSectionIntegration::arrangeFibers(UniaxialMaterial **theMaterials,
					  UniaxialMaterial *theConcrete,
					  UniaxialMaterial *theSteel)
{
  int numFibers = this->getNumFibers();

  int Nfconcrete = Nrings*Nwedges;

  int i;
  for (i = 0; i < Nfconcrete; i++)
    theMaterials[i] = theConcrete;
  for ( ; i < numFibers; i++)
    theMaterials[i] = theSteel;

  return 0;
}

void
RCTunnelSectionIntegration::getFiberLocations(int nFibers, double *yi, double *zi)
{
  static const double pi = 3.141592653589793;
  
  double theta = pi/Nwedges;
  double twoTheta = 2.0*theta;

  int loc = 0;

  //std::ofstream ofs ("fibers.out", std::ofstream::out);


  // 1. Concrete region
  double dr = h/Nrings;
  double rinner = 0.5*d - h;
  double Ainner = rinner*rinner*theta;
  double xinner = 2.0/3.0*rinner*sin(theta)/theta;
  for (int i = 0; i < Nrings; i++) {
    double router = 0.5*d - h + (i+1)*dr;
    double Aouter = router*router*theta;
    double xouter = 2.0/3.0*router*sin(theta)/theta;
    double area = Aouter-Ainner;
    double xbar = (xouter*Aouter-xinner*Ainner)/area;
    double angle = theta;
    for (int j = 0; j < Nwedges; j++) {
      yi[loc] = xbar*cos(angle);
      zi[loc] = xbar*sin(angle);
      //ofs << yi[loc] << ' ' << zi[loc] << endln;
      angle += twoTheta;
      loc++;
    }
    Ainner = Aouter;
    xinner = xouter;
  }

  // 2. Inner steel bars
  double xbar = 0.5*d - h + coverinner;
  theta = pi/Nbarsinner;
  twoTheta = 2.0*theta;
  double angle = theta;
  for (int i = 0; i < Nbarsinner; i++) {
    yi[loc] = xbar*cos(angle);
    zi[loc] = xbar*sin(angle);
    //ofs << yi[loc] << ' ' << zi[loc] << endln;
    angle += twoTheta;
    loc++;
  }

  // 3. Outer steel bars
  xbar = 0.5*d - coverouter;
  theta = pi/Nbarsouter;
  twoTheta = 2.0*theta;
  angle = theta;
  for (int i = 0; i < Nbarsouter; i++) {
    yi[loc] = xbar*cos(angle);
    zi[loc] = xbar*sin(angle);
    //ofs << yi[loc] << ' ' << zi[loc] << endln;
    angle += twoTheta;
    loc++;
  }

  //ofs.close();

  return;
}

void
RCTunnelSectionIntegration::getFiberWeights(int nFibers, double *wt)
{
  static const double pi = 3.141592653589793;

  double theta = pi/Nwedges;
  //double twoTheta = 2*theta;

  int loc = 0;

  //std::ofstream ofs ("areas.out", std::ofstream::out);

  // 1. Concrete region
  double dr = h/Nrings;
  double rinner = 0.5*d - h;
  double Ainner = rinner*rinner*theta;
  for (int i = 0; i < Nrings; i++) {
    double router = 0.5*d - h + (i+1)*dr;
    double Aouter = router*router*theta;
    double area = Aouter-Ainner;
    for (int j = 0; j < Nwedges; j++) {
      wt[loc++] = area;
      //ofs << area << endln;
    }
    Ainner = Aouter;
  }

  // 2. Inner steel bars
  for (int i = 0; i < Nbarsinner; i++) {
    wt[loc++] = Asinner;
    //ofs << Asinner << endln;
  }

  // 3. Outer steel bars
  for (int i = 0; i < Nbarsouter; i++) {
    wt[loc++] = Asouter;
    //ofs << Asouter << endln;
  }

  //ofs.close();

  return;
}

SectionIntegration*
RCTunnelSectionIntegration::getCopy(void)
{
  return new RCTunnelSectionIntegration(d, h, Asinner, Asouter, coverinner, coverouter,
					Nrings, Nwedges, Nbarsinner, Nbarsouter);
}

int
RCTunnelSectionIntegration::setParameter(const char **argv, int argc,
					 Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"d") == 0 || strcmp(argv[0],"D") == 0) {
    param.setValue(d);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"h") == 0) {
    param.setValue(h);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"Asinner") == 0) {
    param.setValue(Asinner);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"Asouter") == 0) {
    param.setValue(Asouter);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"coverinner") == 0) {
    param.setValue(coverinner);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"coverouter") == 0) {
    param.setValue(coverouter);
    return param.addObject(6, this);
  }

  return -1;
}

int
RCTunnelSectionIntegration::updateParameter(int parameterID,
				      Information &info)
{
  switch (parameterID) {
  case 1:
    d = info.theDouble;
    return 0;
  case 2:
    h = info.theDouble;
    return 0;
  case 3:
    Asinner = info.theDouble;
    return 0;
  case 4:
    Asouter = info.theDouble;
    return 0;
  case 5:
    coverinner = info.theDouble;
    return 0;
  case 6:
    coverouter = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
RCTunnelSectionIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
RCTunnelSectionIntegration::getLocationsDeriv(int nFibers, double *dyidh, double *dzidh)
{
  double dddh = 0.0;
  double dhdh = 0.0;
  double dAsinnerdh = 0.0;
  double dAsouterdh = 0.0;
  double dcoverinnerdh = 0.0;
  double dcoverouterdh = 0.0;
  
  if (parameterID == 1) // d
    dddh  = 1.0;
  else if (parameterID == 2) // h
    dhdh = 1.0;
  else if (parameterID == 3) // Asinner
    dAsinnerdh =  1.0;
  else if (parameterID == 4) // Asouter
    dAsouterdh = 1.0;
  else if (parameterID == 5) // coverinner
    dcoverinnerdh =  1.0;
  else if (parameterID == 6) // coverouter
    dcoverouterdh = 1.0;
  else {
    for (int i = 0; i < nFibers; i++) {
      dyidh[i] = 0.0;
      dzidh[i] = 0.0;
    }
    return;
  }

  static const double pi = 3.141592653589793;

  double theta = pi/Nwedges;
  double twoTheta = 2.0*theta;

  int loc = 0;

  // 1. Core region
  double drdh = dhdh/Nrings;
  double rinner = 0.5*dddh - dhdh;
  double Ainner = rinner*rinner*theta;
  double xinner = 2.0/3.0*rinner*sin(theta)/theta;
  for (int i = 0; i < Nrings; i++) {
    double router = 0.5*dddh - dhdh + (i+1)*drdh;
    double Aouter = router*router*theta;
    double xouter = 2.0/3.0*router*sin(theta)/theta;
    double area = Aouter-Ainner;
    double dxbardh = (xouter*Aouter-xinner*Ainner)/area;
    double angle = theta;
    for (int j = 0; j < Nwedges; j++) {
      dyidh[loc] = dxbardh*cos(angle);
      dzidh[loc] = dxbardh*sin(angle);
      angle += twoTheta;
      loc++;
    }
    Ainner = Aouter;
    xinner = xouter;
  }  

  // 2. Inner steel bars
  double dxbardh = 0.5*dddh - dhdh + dcoverinnerdh;
  theta = pi/Nbarsinner;
  twoTheta = 2.0*theta;
  double angle = theta;
  for (int i = 0; i < Nbarsinner; i++) {
    dyidh[loc] = dxbardh*cos(angle);
    dzidh[loc] = dxbardh*sin(angle);
    angle += twoTheta;
    loc++;
  }

  // 3. Outer steel bars
  dxbardh = 0.5*dddh - dcoverouterdh;
  theta = pi/Nbarsouter;
  twoTheta = 2.0*theta;
  angle = theta;
  for (int i = 0; i < Nbarsouter; i++) {
    dyidh[loc] = dxbardh*cos(angle);
    dzidh[loc] = dxbardh*sin(angle);
    angle += twoTheta;
    loc++;
  }

  return;
}

void
RCTunnelSectionIntegration::getWeightsDeriv(int nFibers, double *dwtsdh)
{
  double dddh = 0.0;
  double dhdh = 0.0;
  double dAsinnerdh = 0.0;
  double dAsouterdh = 0.0;
  double dcoverinnerdh = 0.0;
  double dcoverouterdh = 0.0;
  
  if (parameterID == 1) // d
    dddh  = 1.0;
  else if (parameterID == 2) // h
    dhdh = 1.0;
  else if (parameterID == 3) // Asinner
    dAsinnerdh =  1.0;
  else if (parameterID == 4) // Asouter
    dAsouterdh = 1.0;
  else if (parameterID == 5) // coverinner
    dcoverinnerdh =  1.0;
  else if (parameterID == 6) // coverouter
    dcoverouterdh = 1.0;
  else {
    for (int i = 0; i < nFibers; i++)
      dwtsdh[i] = 0.0;
    return;
  }

  static const double pi = 3.141592653589793;

  double theta = pi/Nwedges;

  int loc = 0;

  // 1. Concrete region
  double drdh = dhdh/Nrings;
  double rinner = 0.5*dddh - dhdh;
  double Ainner = rinner*rinner*theta;
  for (int i = 0; i < Nrings; i++) {
    double router = 0.5*dddh - dhdh + (i+1)*drdh;
    double Aouter = router*router*theta;
    double area = Aouter-Ainner;
    for (int j = 0; j < Nwedges; j++)
      dwtsdh[loc++] = area;
    Ainner = Aouter;
  }

  // 2. Inner steel bars
  for (int i = 0; i < Nbarsinner; i++) {
    dwtsdh[loc++] = dAsinnerdh;
  }  

  // 3. Outer steel bars
  for (int i = 0; i < Nbarsouter; i++) {
    dwtsdh[loc++] = dAsouterdh;
  }  

  return;
}

void
RCTunnelSectionIntegration::Print(OPS_Stream &s, int flag)
{
  s << "RC Circular Section" << endln;
  s << " d = " << d;
  s << " h = " << h;
  s << " Asinner = " << Asinner;
  s << " Asouter = " << Asouter;
  s << " coverinner = " << coverinner;
  s << " coverouter = " << coverouter << endln;
  s << " Nrings = " << Nrings;
  s << " Nwedges = " << Nwedges;
  s << " Nbarsinner = " << Nbarsinner;
  s << " Nbarsouter = " << Nbarsouter << endln;

  return;
}

int
RCTunnelSectionIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(10);

  data(0) = d;
  data(1) = h;
  data(2) = Asinner;
  data(3) = Asouter;
  data(4) = coverinner;
  data(5) = coverouter;
  data(6) = Nrings;
  data(7) = Nwedges;
  data(8) = Nbarsinner;
  data(9) = Nbarsouter;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "RCTunnelSectionIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
RCTunnelSectionIntegration::recvSelf(int cTag, Channel &theChannel,
				     FEM_ObjectBroker &theBroker)
{
  static Vector data(10);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "RCTunnelSectionIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  d       = data(0);
  h       = data(1);
  Asinner  = data(2);
  Asouter  = data(3);
  coverinner    = data(4);
  coverinner    = data(5);
  Nrings = (int)data(6);
  Nwedges   = (int)data(7);
  Nbarsinner = (int)data(8);
  Nbarsouter = (int)data(9);

  return 0;
}
