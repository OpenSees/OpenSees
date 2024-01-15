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

#include <RCCircularSectionIntegration.h>
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

void* OPS_RCCircularSection()
{
  if (OPS_GetNumRemainingInputArgs() < 11) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: section RCCircularSection tag? coreTag? coverTag? steelTag? d? cover? As? NringsCore? NringsCover? Nwedges? Nsteel? <-GJ GJ?> <or> <-torsion matTag?>\n";
    return 0;
  }

  int ndm = OPS_GetNDM();
  
  int idata[8];
  double ddata[3];
  
  int numdata = 4;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid section RCCircularSection input\n";
    return 0;
  }
  
  numdata = 3;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "WARNING invalid section RCCircularSection input\n";
    return 0;
  }
  
  numdata = 4;
  if (OPS_GetIntInput(&numdata, &idata[4]) < 0) {
    opserr << "WARNING invalid section RCCircularSection input\n";
    return 0;
  }
  
  int tag=idata[0], coreTag=idata[1], coverTag=idata[2], steelTag=idata[3];
  double d=ddata[0], cover=ddata[1], As=ddata[2];
  int ncore=idata[4], ncover=idata[5], nwedge=idata[6], nsteel=idata[7];
  
  UniaxialMaterial *theCore = OPS_getUniaxialMaterial(coreTag);
  
  if (theCore == 0) {
    opserr << "WARNING uniaxial material does not exist\n";
    opserr << "material: " << coreTag; 
    opserr << "\nRCCircularSection section: " << tag << endln;
    return 0;
  }
  
  UniaxialMaterial *theCover = OPS_getUniaxialMaterial(coverTag);
  
  if (theCover == 0) {
    opserr << "WARNING uniaxial material does not exist\4n";
    opserr << "material: " << coverTag; 
    opserr << "\nRCCircularSection section: " << tag << endln;
    return 0;
  }
  
  UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(steelTag);
  
  if (theSteel == 0) {
    opserr << "WARNING uniaxial material does not exist\n";
    opserr << "material: " << steelTag; 
    opserr << "\nRCCircularSection section: " << tag << endln;
    return 0;
  }
  
  RCCircularSectionIntegration rcsect(d, As, cover, ncore, ncover, nwedge, nsteel);
  
  int numFibers = rcsect.getNumFibers();
  
  UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];
  
  rcsect.arrangeFibers(theMats, theCore, theCover, theSteel);
  
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
    opserr << "WARNING torsion not specified for RCCircularSection\n";
    opserr << "\nRCCircularSection section: " << tag << endln;
    return 0;
  }
  
  // Parsing was successful, allocate the section
  SectionForceDeformation* theSection = 0;
  if (ndm == 2)
    theSection = new FiberSection2d(tag, numFibers, theMats, rcsect);
  if (ndm == 3)
    theSection = new FiberSection3d(tag, numFibers, theMats, rcsect, *torsion);
  
  if (deleteTorsion)
    delete torsion;
  
  delete [] theMats;
  
  return theSection;
}

RCCircularSectionIntegration::RCCircularSectionIntegration(double D,
					   double AS,
					   double COVER,
					   int NCORE,
					   int NCOVER,
					   int NWEDGES,
					   int NFS):
  SectionIntegration(SECTION_INTEGRATION_TAG_RCCIRCULAR),
  d(D), As(AS), cover(COVER),
  NringsCore(NCORE), NringsCover(NCOVER), Nwedges(NWEDGES), Nsteel(NFS), 
  parameterID(0)
{
  /*
  if (NringsCore < 1)
    NringsCore = 1;

  if (NringsCover < 1)
    NringsCover = 1;

  if (Nwedges < 2)
    Nwedges = 2;

  if (Nsteel < 1)
    Nsteel = 1;
  */
}

RCCircularSectionIntegration::RCCircularSectionIntegration():
  SectionIntegration(SECTION_INTEGRATION_TAG_RCCIRCULAR),
  d(0.0), As(0.0), cover(0.0),
  NringsCore(1), NringsCover(1), Nwedges(2), Nsteel(1),
  parameterID(0)
{
  
}

RCCircularSectionIntegration::~RCCircularSectionIntegration()
{
  
}

int
RCCircularSectionIntegration::getNumFibers(FiberType type)
{
  if (type == steel)
    return Nsteel;
  if (type == concrete)
    return (NringsCore+NringsCover)*Nwedges;
  if (type == all)
    return (NringsCore+NringsCover)*Nwedges + Nsteel;

  return 0;
}

int
RCCircularSectionIntegration::arrangeFibers(UniaxialMaterial **theMaterials,
					    UniaxialMaterial *theCore,
					    UniaxialMaterial *theCover,
					    UniaxialMaterial *theSteel)
{
  int numFibers = this->getNumFibers();

  int Nfcore = NringsCore*Nwedges;

  int i;
  for (i = 0; i < Nfcore; i++)
    theMaterials[i] = theCore;
  for ( ; i < numFibers-Nsteel; i++)
    theMaterials[i] = theCover;
  for ( ; i < numFibers; i++)
    theMaterials[i] = theSteel;

  return 0;
}

void
RCCircularSectionIntegration::getFiberLocations(int nFibers, double *yi, double *zi)
{
  static const double pi = 3.141592653589793;
  
  double theta = pi/Nwedges;
  double twoTheta = 2.0*theta;

  int loc = 0;

  //std::ofstream ofs ("fibers.out", std::ofstream::out);


  // 1. Core region
  double dr = (0.5*d-cover)/NringsCore;
  double rinner = 0.0;
  double Ainner = 0.0;
  double xinner = 0.0;
  for (int i = 0; i < NringsCore; i++) {
    double router = (i+1)*dr;
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

  // 2. Cover region
  dr = cover/NringsCover;
  rinner = 0.5*d - cover;
  // use xinner from above
  Ainner = rinner*rinner*theta;
  for (int i = 0; i < NringsCover; i++) {
    double router = 0.5*d - cover + (i+1)*dr;
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

  // 3. Steel bars
  double xbar = 0.5*d - cover;
  theta = pi/Nsteel;
  twoTheta = 2.0*theta;
  double angle = theta;
  for (int i = 0; i < Nsteel; i++) {
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
RCCircularSectionIntegration::getFiberWeights(int nFibers, double *wt)
{
  static const double pi = 3.141592653589793;

  double theta = pi/Nwedges;
  //double twoTheta = 2*theta;

  int loc = 0;

  //std::ofstream ofs ("areasCircular.out", std::ofstream::out);

  // 1. Core region
  double dr = (0.5*d-cover)/NringsCore;
  double rinner = 0.0;
  double Ainner = 0.0;
  for (int i = 0; i < NringsCore; i++) {
    double router = (i+1)*dr;
    double Aouter = router*router*theta;
    double area = Aouter-Ainner;
    for (int j = 0; j < Nwedges; j++) {
      wt[loc++] = area;
      //ofs << area << endln;
    }
    Ainner = Aouter;
  }

  // 2. Cover region
  dr = cover/NringsCover;
  rinner = 0.5*d - cover;
  Ainner = rinner*rinner*theta;
  for (int i = 0; i < NringsCover; i++) {
    double router = 0.5*d - cover + (i+1)*dr;
    double Aouter = router*router*theta;
    double area = Aouter-Ainner;
    for (int j = 0; j < Nwedges; j++) {
      wt[loc++] = area;
      //ofs << area << endln;
    }
    Ainner = Aouter;
  }

  // 3. Steel bars
  for (int i = 0; i < Nsteel; i++) {
    wt[loc++] = As;
    //ofs << As << endln;
  }

  //ofs.close();

  return;
}

SectionIntegration*
RCCircularSectionIntegration::getCopy(void)
{
  return new RCCircularSectionIntegration(d, As, cover, NringsCore,
					  NringsCover, Nwedges, Nsteel);
}

int
RCCircularSectionIntegration::setParameter(const char **argv, int argc,
				   Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"d") == 0 || strcmp(argv[0],"D") == 0) {
    param.setValue(d);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"As") == 0) {
    param.setValue(As);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"cover") == 0) {
    param.setValue(cover);
    return param.addObject(6, this);
  }

  return -1;
}

int
RCCircularSectionIntegration::updateParameter(int parameterID,
				      Information &info)
{
  switch (parameterID) {
  case 1:
    d = info.theDouble;
    return 0;
  case 5:
    As = info.theDouble;
    return 0;
  case 6:
    cover = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
RCCircularSectionIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
RCCircularSectionIntegration::getLocationsDeriv(int nFibers, double *dyidh, double *dzidh)
{

  double dddh  = 0.0;
  double dAsdh = 0.0;
  double dcoverdh = 0.0;
  
  if (parameterID == 1) // d
    dddh  = 1.0;
  else if (parameterID == 5) // As
    dAsdh = 1.0;
  else if (parameterID == 6) // cover
    dcoverdh =  1.0;
  else {
    for (int i = 0; i < nFibers; i++) {
      dyidh[i] = 0.0;
      if (dzidh != 0)
	dzidh[i] = 0.0;
    }
    return;
  }

  static const double pi = 3.141592653589793;

  double theta = pi/Nwedges;
  double twoTheta = 2.0*theta;

  int loc = 0;

  // 1. Core region
  double dr = (0.5*d-cover)/NringsCore;
  double ddrdh = (0.5*dddh-dcoverdh)/NringsCore;
  double rinner = 0.0;
  double drinnerdh = 0.0;
  double Ainner = 0.0;
  double dAinnerdh = 0.0;
  double xinner = 0.0;
  double dxinnerdh = 0.0;
  for (int i = 0; i < NringsCore; i++) {
    double router = (i+1)*dr;
    double drouterdh = (i+1)*ddrdh;
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
  
  // 2. Cover region
  dr = cover/NringsCover;
  ddrdh = dcoverdh/NringsCover;
  rinner = 0.5*d - cover;
  drinnerdh = 0.5*dddh - dcoverdh;
  // use xinner and dxinnerdh from above
  Ainner = rinner*rinner*theta;
  dAinnerdh = 2*rinner*drinnerdh*theta;
  for (int i = 0; i < NringsCover; i++) {
    double router = 0.5*d - cover + (i+1)*dr;
    double drouterdh = 0.5*dddh - dcoverdh + (i+1)*ddrdh;
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

  // 3. Steel bars
  double dxbardh = 0.5*dddh-dcoverdh;
  theta = pi/Nsteel;
  twoTheta = 2.0*theta;
  double angle = theta;
  for (int i = 0; i < Nsteel; i++) {
    dyidh[loc] = dxbardh*cos(angle);
    if (dzidh != 0)
      dzidh[loc] = dxbardh*sin(angle);
    angle += twoTheta;
    loc++;
  }

  return;
}

void
RCCircularSectionIntegration::getWeightsDeriv(int nFibers, double *dwtsdh)
{
  double dddh  = 0.0;
  double dAsdh = 0.0;
  double dcoverdh = 0.0;
  
  if (parameterID == 1) { // d 
    dddh  = 1.0;
  }
  else if (parameterID == 5) // As
    dAsdh = 1.0;
  else if (parameterID == 6) // cover
    dcoverdh =  1.0;
  else {
    for (int i = 0; i < nFibers; i++)
      dwtsdh[i] = 0.0;
    return;
  }

  static const double pi = 3.141592653589793;

  double theta = pi/Nwedges;

  int loc = 0;

  // 1. Core region
  double dr = (0.5*d-cover)/NringsCore;
  double ddrdh = (0.5*dddh-dcoverdh)/NringsCore;
  double rinner = 0.0;
  double drinnerdh = 0.0;
  double Ainner = 0.0;
  double dAinnerdh = 0.0;
  for (int i = 0; i < NringsCore; i++) {
    double router = (i+1)*dr;
    double drouterdh = (i+1)*ddrdh;
    double Aouter = router*router*theta;
    double dAouterdh = 2*router*drouterdh*theta;
    double area = Aouter-Ainner;
    double dareadh = dAouterdh-dAinnerdh;
    for (int j = 0; j < Nwedges; j++)
      dwtsdh[loc++] = dareadh;
    Ainner = Aouter;
    dAinnerdh = dAouterdh;
  }

  // 2. Cover region
  dr = cover/NringsCover;
  ddrdh = dcoverdh/NringsCover;
  rinner = 0.5*d - cover;
  drinnerdh = 0.5*dddh - dcoverdh;
  Ainner = rinner*rinner*theta;
  dAinnerdh = 2*rinner*drinnerdh*theta;
  for (int i = 0; i < NringsCover; i++) {
    double router = 0.5*d - cover + (i+1)*dr;
    double drouterdh = 0.5*dddh - dcoverdh + (i+1)*ddrdh;
    double Aouter = router*router*theta;
    double dAouterdh = 2*router*drouterdh*theta;
    double area = Aouter-Ainner;
    double dareadh = dAouterdh-dAinnerdh;
    for (int j = 0; j < Nwedges; j++)
      dwtsdh[loc++] = dareadh;
    Ainner = Aouter;
    dAinnerdh = dAouterdh;
  }

  // 3. Steel bars
  for (int i = 0; i < Nsteel; i++) {
    dwtsdh[loc++] = dAsdh;
  }  

  return;
}

void
RCCircularSectionIntegration::Print(OPS_Stream &s, int flag)
{
  s << "RC Circular Section" << endln;
  s << " d = "  << d;
  s << " As = " << As;
  s << " cover = " << cover << endln;
  s << " NringsCore = " << NringsCore;
  s << " NringsCover = " << NringsCover;
  s << " Nwedges = " << Nwedges;
  s << " Nsteel = " << Nsteel << endln;

  return;
}

int
RCCircularSectionIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(7);

  data(0) = d;
  data(1) = As;
  data(2) = cover;
  data(3) = NringsCore;
  data(4) = NringsCover;
  data(5) = Nwedges;
  data(6) = Nsteel;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "RCCircularSectionIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
RCCircularSectionIntegration::recvSelf(int cTag, Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
{
  static Vector data(7);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "RCCircularSectionIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  d       = data(0);
  As       = data(1);
  cover    = data(2);
  NringsCore = (int)data(3);
  NringsCover   = (int)data(4);
  Nwedges   = (int)data(5);
  Nsteel  = (int)data(6);

  return 0;
}
