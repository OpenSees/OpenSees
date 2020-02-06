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

#include <elementAPI.h>
#include <UniaxialMaterial.h>
#include <ElasticMaterial.h>
#include <NDMaterial.h>
#include <FiberSection3d.h>
#include <NDFiberSection3d.h>

void* OPS_TubeSection()
{
  if (OPS_GetNumRemainingInputArgs() < 6) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: section Tube tag? matTag? D? t? nfw? nfr? <-nd shape?>" << endln;
    return 0;
  }
  
  int tag, matTag;
  double D, t;
  int nfw, nfr;
  
  SectionForceDeformation* theSection = 0;
  
  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid section Tube tag" << endln;
    return 0;
  }
  
  if (OPS_GetIntInput(&numdata, &matTag) < 0) {
    opserr << "WARNING invalid section Tube matTag" << endln;
    return 0;
  }
  
  if (OPS_GetDoubleInput(&numdata, &D) < 0) {
    opserr << "WARNING invalid D" << endln;
    opserr << "Tube section: " << tag << endln;
    return 0;
  }
  
  if (OPS_GetDoubleInput(&numdata, &t) < 0) {
    opserr << "WARNING invalid t" << endln;
    opserr << "Tube section: " << tag << endln;
    return 0;
  }
  
  if (OPS_GetIntInput(&numdata, &nfw) < 0) {
    opserr << "WARNING invalid nfw" << endln;
    opserr << "Tube section: " << tag << endln;
    return 0;
  }
  
  if (OPS_GetIntInput(&numdata, &nfr) < 0) {
    opserr << "WARNING invalid nfr" << endln;
    opserr << "Tube section: " << tag << endln;
    return 0;
  }
  
  TubeSectionIntegration tubesect(D, t, nfw, nfr);
  
  int numFibers = tubesect.getNumFibers();
  
  if (OPS_GetNumRemainingInputArgs() > 0) {
    
    double shape = 1.0;
    if (OPS_GetNumRemainingInputArgs() > 1) {
      if (OPS_GetDoubleInput(&numdata, &shape) < 0) {
	opserr << "WARNING invalid shape" << endln;
	opserr << "Tube section: " << tag << endln;
	return 0;
      }
    }
    
    NDMaterial *theSteel = OPS_getNDMaterial(matTag);
    
    if (theSteel == 0) {
      opserr << "WARNING ND material does not exist\n";
      opserr << "material: " << matTag;
      opserr << "\nTube section: " << tag << endln;
      return 0;
    }
    
    NDMaterial **theMats = new NDMaterial *[numFibers];
    
    tubesect.arrangeFibers(theMats, theSteel);
    
    // Parsing was successful, allocate the section
    theSection = 0;
    if (OPS_GetNumRemainingInputArgs() > 0) {
      const char* flag = OPS_GetString();
      if (strcmp(flag,"-nd") == 0) {
	theSection = new NDFiberSection3d(tag, numFibers, theMats, tubesect, shape);
      } else if (strcmp(flag,"-ndWarping") == 0) {
	opserr << "TubeSection -- not implemented yet for fiber warping section" << endln;
	theSection = 0;
	//theSection = new NDFiberSectionWarping2d(tag, numFibers, theMats, tubesect, shape);
      }
    }
    delete [] theMats;
  }
  else {
    UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(matTag);
    
    if (theSteel == 0) {
      opserr << "WARNING uniaxial material does not exist\n";
      opserr << "material: " << matTag;
      opserr << "\nTube section: " << tag << endln;
      return 0;
    }
    
    UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];
    
    tubesect.arrangeFibers(theMats, theSteel);


    UniaxialMaterial *torsion = 0;
    if (OPS_GetNumRemainingInputArgs() < 2) {
      opserr << "WARNING torsion not specified for TubeSection\n";
      opserr << "Use either -GJ $GJ or -torsion $matTag\n";
      opserr << "\nTubeSection: " << tag << endln;
      return 0;
    }
    const char* opt = OPS_GetString();
    numdata = 1;
    bool deleteTorsion = false;
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
    if (torsion == 0) {
      opserr << "WARNING torsion not speified for TubeSection\n";
      opserr << "\nTubeSection section: " << tag << endln;
      return 0;
    }
    
    // Parsing was successful, allocate the section
    theSection = new FiberSection3d(tag, numFibers, theMats, tubesect, *torsion);
    
    delete [] theMats;
  }
  
  return theSection;
}



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
