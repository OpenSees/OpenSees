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

#include <elementAPI.h>
#include <UniaxialMaterial.h>
#include <ElasticMaterial.h>
#include <FiberSection3d.h>
#include <FiberSection2d.h>

void* OPS_RCTunnelSection()
{
  if (OPS_GetNumRemainingInputArgs() < 13) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: section RCTunnelSection tag? concreteTag? steelTag? d? h? coverinner? coverouter? Asinner? Asouter? Nrings? Nwedges? Nbarsinner? Nbarsouter? <-GJ GJ?> <or> <-torsion matTag?>\n";
    return 0;
  }

  int ndm = OPS_GetNDM();
  
  int idata[8];
  double ddata[6];
  
  int numdata = 3;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid section RCTunnelSection input\n";
    return 0;
  }
  
  numdata = 6;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "WARNING invalid section RCTunnelSection input\n";
    return 0;
  }
  
  numdata = 4;
  if (OPS_GetIntInput(&numdata, &idata[4]) < 0) {
    opserr << "WARNING invalid section RCTunnelSection input\n";
    return 0;
  }
  
  int tag=idata[0], concreteTag=idata[1], steelTag=idata[2];
  double d=ddata[0], h=ddata[1], coverinner=ddata[2], coverouter=ddata[3], Asinner=ddata[4], Asouter=ddata[5];
  int nring=idata[3], nwedge=idata[4], nbarsinner=idata[5], nbarsouter=idata[6];
  
  UniaxialMaterial *theConcrete = OPS_getUniaxialMaterial(concreteTag);
  
  if (theConcrete == 0) {
    opserr << "WARNING uniaxial material does not exist\n";
    opserr << "material: " << concreteTag; 
    opserr << "\nRCTunnelSection section: " << tag << endln;
    return 0;
  }
  
  UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(steelTag);
  
  if (theSteel == 0) {
    opserr << "WARNING uniaxial material does not exist\n";
    opserr << "material: " << steelTag; 
    opserr << "\nRCTunnelSection section: " << tag << endln;
    return 0;
  }
  
  RCTunnelSectionIntegration rcsect(d, h, Asinner, Asouter, coverinner, coverouter,
				    nring, nwedge, nbarsinner, nbarsouter);
  
  int numFibers = rcsect.getNumFibers();
  
  UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];
  
  rcsect.arrangeFibers(theMats, theConcrete, theSteel);
  
  UniaxialMaterial *torsion = 0;
  numdata = 1;
  bool deleteTorsion = false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    if (strcmp(opt, "-GJ") == 0) {
      double GJ;
      if (OPS_GetDoubleInput(&numdata, &GJ) < 0) {
	opserr << "WARNING: failed to read GJ\n";
	return 0;
      }
      torsion = new ElasticMaterial(0,GJ);
      deleteTorsion = true;
    }
    if (strcmp(opt, "-torsion") == 0) {
      int torsionTag;
      if (OPS_GetIntInput(&numdata, &torsionTag) < 0) {
	opserr << "WARNING: failed to read torsion\n";
	return 0;
      }
      torsion = OPS_getUniaxialMaterial(torsionTag);
    }
  }
  if (ndm == 3 && torsion == 0) {
    opserr << "WARNING torsion not speified for RCCircularSection\n";
    opserr << "\nRCTunnelSection section: " << tag << endln;
    return 0;
  }
  
  // Parsing was successful, allocate the section
  SectionForceDeformation* theSection = 0;
  if (ndm == 2)
    theSection = new FiberSection2d(tag, numFibers, theMats, rcsect);
  if (ndm == 3)
    theSection = new FiberSection3d(tag, numFibers, theMats, rcsect, *torsion);  
  
  delete [] theMats;
  if (deleteTorsion)
    delete torsion;
  
  return theSection;
}

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
  double rinner = 0.5*d;
  double Ainner = rinner*rinner*theta;
  double xinner = 2.0/3.0*rinner*sin(theta)/theta;
  for (int i = 0; i < Nrings; i++) {
    double router = 0.5*d + (i+1)*dr;
    double Aouter = router*router*theta;
    double xouter = 2.0/3.0*router*sin(theta)/theta;
    double area = Aouter-Ainner;
    double xbar = (xouter*Aouter-xinner*Ainner)/area;
    double angle = theta;
    for (int j = 0; j < Nwedges; j++) {
      yi[loc] = xbar*cos(angle);
      if (zi != 0)
	zi[loc] = xbar*sin(angle);
      //ofs << yi[loc] << ' ' << zi[loc] << endln;
      angle += twoTheta;
      loc++;
    }
    Ainner = Aouter;
    xinner = xouter;
  }

  // 2. Inner steel bars
  double xbar = 0.5*d + coverinner;
  theta = pi/Nbarsinner;
  twoTheta = 2.0*theta;
  double angle = theta;
  for (int i = 0; i < Nbarsinner; i++) {
    yi[loc] = xbar*cos(angle);
    if (zi != 0)
      zi[loc] = xbar*sin(angle);
    //ofs << yi[loc] << ' ' << zi[loc] << endln;
    angle += twoTheta;
    loc++;
  }

  // 3. Outer steel bars
  xbar = 0.5*d + h - coverouter;
  theta = pi/Nbarsouter;
  twoTheta = 2.0*theta;
  angle = theta;
  for (int i = 0; i < Nbarsouter; i++) {
    yi[loc] = xbar*cos(angle);
    if (zi != 0)
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

  //std::ofstream ofs ("areasTunnel.out", std::ofstream::out);

  // 1. Concrete region
  double dr = h/Nrings;
  double rinner = 0.5*d;
  double Ainner = rinner*rinner*theta;
  for (int i = 0; i < Nrings; i++) {
    double router = 0.5*d + (i+1)*dr;
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
    for (int i = 0; i < nFibers; i++)
      dyidh[i] = 0.0;
    if (dzidh != 0) {
      for (int i = 0; i < nFibers; i++)      
	dzidh[i] = 0.0;
    }
    return;
  }

  static const double pi = 3.141592653589793;

  // theta, Nrings, and Nwedges are constant
  // when taking derivatives
  double theta = pi/Nwedges;
  double twoTheta = 2.0*theta;

  int loc = 0;

  // 1. Concrete region
  double dr = h/Nrings;
  double ddrdh = dhdh/Nrings;
  double rinner = 0.5*d;
  double drinnerdh = 0.5*dddh;
  double Ainner = rinner*rinner*theta;
  double dAinnerdh = 2*rinner*drinnerdh*theta;
  double xinner = 2.0/3.0*rinner*sin(theta)/theta;
  double dxinnerdh = 2.0/3.0*drinnerdh*sin(theta)/theta;
  for (int i = 0; i < Nrings; i++) {
    double router = 0.5*d + (i+1)*dr;
    double drouterdh = 0.5*dddh + (i+1)*ddrdh;
    double Aouter = router*router*theta;
    double dAouterdh = 2*router*drouterdh*theta;
    double xouter = 2.0/3.0*router*sin(theta)/theta;
    double dxouterdh = 2.0/3.0*drouterdh*sin(theta)/theta;
    double area = Aouter-Ainner;
    double dareadh = dAouterdh-dAinnerdh;
    double xbar = (xouter*Aouter-xinner*Ainner)/area;
    double dxbardh = (area*(xouter*dAouterdh+dxouterdh*Aouter-xinner*dAinnerdh-dxinnerdh*Ainner)-(xouter*Aouter-xinner*Ainner)*dareadh)/(area*area);
    double angle = theta;
    for (int j = 0; j < Nwedges; j++) {
      dyidh[loc] = dxbardh*cos(angle);
      if (dzidh != 0)
	dzidh[loc] = dxbardh*sin(angle);
      angle += twoTheta;
      loc++;
    }
    Ainner = Aouter;
    dAinnerdh = dAouterdh;
    xinner = xouter;
    dxinnerdh = dxouterdh;
  }  

  // 2. Inner steel bars
  double dxbardh = 0.5*dddh + dcoverinnerdh;
  theta = pi/Nbarsinner;
  twoTheta = 2.0*theta;
  double angle = theta;
  for (int i = 0; i < Nbarsinner; i++) {
    dyidh[loc] = dxbardh*cos(angle);
    if (dzidh != 0)
      dzidh[loc] = dxbardh*sin(angle);
    angle += twoTheta;
    loc++;
  }

  // 3. Outer steel bars
  dxbardh = 0.5*dddh + dhdh - dcoverouterdh;
  theta = pi/Nbarsouter;
  twoTheta = 2.0*theta;
  angle = theta;
  for (int i = 0; i < Nbarsouter; i++) {
    dyidh[loc] = dxbardh*cos(angle);
    if (dzidh != 0)
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

  // theta, Nrings, and Nwedges are constant
  // when taking derivatives
  double theta = pi/Nwedges;

  int loc = 0;

  // 1. Concrete region
  double dr = h/Nrings;
  double ddrdh = dhdh/Nrings;
  double rinner = 0.5*d;
  double drinnerdh = 0.5*dddh;
  double Ainner = rinner*rinner*theta;
  double dAinnerdh = 2*rinner*drinnerdh*theta;
  for (int i = 0; i < Nrings; i++) {
    double router = 0.5*d + (i+1)*dr;
    double drouterdh = 0.5*dddh + (i+1)*ddrdh;
    double Aouter = router*router*theta;
    double dAouterdh = 2*router*drouterdh*theta;
    double area = Aouter-Ainner;
    double dareadh = dAouterdh-dAinnerdh;
    for (int j = 0; j < Nwedges; j++)
      dwtsdh[loc++] = dareadh;
    Ainner = Aouter;
    dAinnerdh = dAouterdh;
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
  s << " D_inner = " << d;
  s << " h = " << h;
  s << " As_inner = " << Asinner;
  s << " As_outer = " << Asouter;
  s << " cover_inner = " << coverinner;
  s << " cover_outer = " << coverouter << endln;
  s << " Nrings = " << Nrings;
  s << " Nwedges = " << Nwedges;
  s << " Nbars_inner = " << Nbarsinner;
  s << " Nbars_outer = " << Nbarsouter << endln;

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
  coverouter    = data(5);
  Nrings = (int)data(6);
  Nwedges   = (int)data(7);
  Nbarsinner = (int)data(8);
  Nbarsouter = (int)data(9);

  return 0;
}
