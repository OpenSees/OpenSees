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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/integration/HSSSectionIntegration.cpp,v $

#include <HSSSectionIntegration.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

#include <elementAPI.h>
#include <UniaxialMaterial.h>
#include <ElasticMaterial.h>
#include <NDMaterial.h>
#include <FiberSection2d.h>
#include <FiberSection3d.h>
#include <NDFiberSection2d.h>
#include <NDFiberSection3d.h>
#include <NDFiberSectionWarping2d.h>

void* OPS_HSSSection()
{
  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: section HSS tag? matTag? h? b? t? nfh? nfb? <-nd> <-shape shape?> <-GJ GJ?> <-torsion tag?>" << endln;
    return 0;
  }
  
  int ndm = OPS_GetNDM();
  
  int tag, matTag;
  double h, b, t;
  int nfh, nfb;
  int nft = 1;
  
  SectionForceDeformation* theSection = 0;
  
  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid section HSS tag" << endln;
    return 0;
  }
  
  if (OPS_GetIntInput(&numdata, &matTag) < 0) {
    opserr << "WARNING invalid section HSS matTag" << endln;
    return 0;
  }
  
  if (OPS_GetDoubleInput(&numdata, &h) < 0) {
    opserr << "WARNING invalid h" << endln;
    opserr << "HSS section: " << tag << endln;
    return 0;
  }
  
  if (OPS_GetDoubleInput(&numdata, &b) < 0) {
    opserr << "WARNING invalid b" << endln;
    opserr << "HSS section: " << tag << endln;
    return 0;
  }
  
  if (OPS_GetDoubleInput(&numdata, &t) < 0) {
    opserr << "WARNING invalid t" << endln;
    opserr << "HSS section: " << tag << endln;
    return 0;
  }
  
  if (OPS_GetIntInput(&numdata, &nfh) < 0) {
    opserr << "WARNING invalid nfh" << endln;
    opserr << "HSS section: " << tag << endln;
    return 0;
  }
  
  if (OPS_GetIntInput(&numdata, &nfb) < 0) {
    opserr << "WARNING invalid nfb" << endln;
    opserr << "HSS section: " << tag << endln;
    return 0;
  }
  
  HSSSectionIntegration hsssect(h, b, t, nfh, nfb, nft);
  
  int numFibers = hsssect.getNumFibers();
  bool isND = false; double shape = 1.0;
  UniaxialMaterial *torsion = 0;
  bool deleteTorsion = false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();
    // read <-nd>
    if (strcmp(flag,"-nd") == 0)
      isND = true;
    // read <-shape shape>
    if (strcmp(flag,"-shape") == 0 && OPS_GetNumRemainingInputArgs() > 0) {
      if (OPS_GetDoubleInput(&numdata, &shape) < 0) {
	opserr << "WARNING invalid shape" << endln;
	opserr << "HSS section: " << tag << endln;
	return 0;
      }
      isND = true;
    }
    // read <-GJ GJ>
    if (strcmp(flag,"-GJ") == 0 && OPS_GetNumRemainingInputArgs() > 0) {
      double GJ;
      if (OPS_GetDoubleInput(&numdata, &GJ) < 0) {
	opserr << "WARNING: failed to read GJ\n";
	return 0;
      }
      torsion = new ElasticMaterial(0,GJ);
      deleteTorsion = true;	    
    }
    // read <-torsion tag>
    if (strcmp(flag,"-torsion") == 0 && OPS_GetNumRemainingInputArgs() > 0) {
      int torsionTag;
      if (OPS_GetIntInput(&numdata, &torsionTag) < 0) {
	opserr << "WARNING: failed to read torsion\n";
	return 0;
      }
      torsion = OPS_getUniaxialMaterial(torsionTag);	    
    }
  }
  
  if (isND) {
    NDMaterial *theSteel = OPS_getNDMaterial(matTag);
    if (theSteel == 0) {
      opserr << "WARNING ND material does not exist\n";
      opserr << "material: " << matTag;
      opserr << "\nHSS section: " << tag << endln;
      return 0;
    }
    
    NDMaterial **theMats = new NDMaterial *[numFibers];
    hsssect.arrangeFibers(theMats, theSteel);

    if (ndm == 2)
      theSection = new NDFiberSection2d(tag, numFibers, theMats, hsssect);
    if (ndm == 3)
      theSection = new NDFiberSection3d(tag, numFibers, theMats, hsssect, shape);
    
    delete [] theMats;
  }

  else {
    UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(matTag);
    if (theSteel == 0) {
      opserr << "WARNING uniaxial material does not exist\n";
      opserr << "material: " << matTag;
      opserr << "\nHSS section: " << tag << endln;
      return 0;
    }
    
    if (ndm == 3 && torsion == 0) {
      opserr << "WARNING torsion not specified for FiberSection\n";
      opserr << "\nHSS section: " << tag << endln;
      return 0;
    }
    
    UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];
    hsssect.arrangeFibers(theMats, theSteel);
    
    if (ndm == 2)
      theSection = new FiberSection2d(tag, numFibers, theMats, hsssect);
    if (ndm == 3)
      theSection = new FiberSection3d(tag, numFibers, theMats, hsssect, *torsion);
    
    if (deleteTorsion)
      delete torsion;
    
    delete [] theMats;
  }
  
  return theSection;  
}  

HSSSectionIntegration::HSSSectionIntegration(double H, double B, double T,
					     int NFH, int NFB, int NFT):
  SectionIntegration(SECTION_INTEGRATION_TAG_HSS),
  h(H), b(B), t(T), Nfh(NFH), Nfb(NFB), Nft(1), parameterID(0)
{
  
}

HSSSectionIntegration::HSSSectionIntegration():
  SectionIntegration(SECTION_INTEGRATION_TAG_HSS),
  h(0.0), b(0.0), t(0.0), Nfh(0), Nfb(0), Nft(0), parameterID(0)
{
  
}

HSSSectionIntegration::~HSSSectionIntegration()
{
  
}

int
HSSSectionIntegration::getNumFibers(FiberType type)
{
  return 2*Nfh + 2*Nfb + 4*Nft*Nft;
}

int
HSSSectionIntegration::arrangeFibers(UniaxialMaterial **theMaterials,
				     UniaxialMaterial *theSteel)
{
  int numFibers = this->getNumFibers();

  for (int i = 0; i < numFibers; i++)
    theMaterials[i] = theSteel;
  
  return 0;
}

int
HSSSectionIntegration::arrangeFibers(NDMaterial **theMaterials,
				     NDMaterial *theSteel)
{
  int numFibers = this->getNumFibers();

  for (int i = 0; i < numFibers; i++)
    theMaterials[i] = theSteel;
  
  return 0;
}

void
HSSSectionIntegration::getFiberLocations(int nFibers, double *yi, double *zi)
{
  double hw = h - 2*t;
  double bw = b - 2*t;  
  
  int i, loc;
  
  double yIncr  = hw/Nfh;
  double yStart = 0.5 * (hw-yIncr);
  
  for (loc = 0, i = 0; loc < Nfh; loc++, i++) {
    yi[loc]     = yStart - yIncr*i;
    yi[loc+Nfh] = yi[loc];
  }
  if (zi != 0) {
    for (loc = 0; loc < Nfh; loc++) {
      zi[loc]     = 0.5*(bw + t);
      zi[loc+Nfh] = -zi[loc];
    }
  }
  
  for (loc = 2*Nfh, i = 0; i < Nfb; loc++, i++) {
    yi[loc]     = 0.5*(hw + t);
    yi[loc+Nfb] = -yi[loc];
  }
  if (zi != 0) {
    double zIncr  = bw/Nfb;
    double zStart = 0.5 * (bw-zIncr);
    for (loc = 2*Nfh, i = 0; i < Nfb; loc++, i++) {
      zi[loc]     = zStart - zIncr*i;
      zi[loc+Nfb] = zi[loc];
    }
  }

  // Four corners
  loc = 2*Nfh + 2*Nfb;
  double yy = 0.5*(hw+t);
  yi[loc++] =  yy;
  yi[loc++] = -yy;
  yi[loc++] = -yy;
  yi[loc++] =  yy;
  if (zi != 0) {
    loc = 2*Nfh + 2*Nfb;
    double zz = 0.5*(bw+t);
    zi[loc++] =  zz;
    zi[loc++] =  zz;
    zi[loc++] = -zz;
    zi[loc++] = -zz;
  }

  return;
}

void
HSSSectionIntegration::getFiberWeights(int nFibers, double *wt)
{
  double hw = h - 2*t;
  double bw = b - 2*t;

  // Assuming Nft = 1
  double a_h = hw*t/(Nfh);
  double a_b = bw*t/(Nfb);
  
  int i, loc;
  
  for (loc = 0; loc < Nfh; loc++) {
    wt[loc]     = a_h;
    wt[loc+Nfh] = a_h;
  }

  for (loc = 2*Nfh, i = 0; i < Nfb; i++, loc++) {
    wt[loc]     = a_b;
    wt[loc+Nfb] = a_b;    
  }

  // Four corners
  double a_t = t*t;
  loc = 2*Nfh + 2*Nfb;
  for (i = 0; i < 4; i++, loc++) {
    wt[loc] = a_t;
  }
  
  return;
}

SectionIntegration*
HSSSectionIntegration::getCopy(void)
{
  HSSSectionIntegration* theCopy = 
    new HSSSectionIntegration(h, b, t, Nfh, Nfb, Nft);

  theCopy->parameterID = parameterID;

  return theCopy;
}

int
HSSSectionIntegration::setParameter(const char **argv, int argc,
					   Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"h") == 0) {
    param.setValue(h);
    return param.addObject(1, this);    
  }
  if (strcmp(argv[0],"b") == 0) {
    param.setValue(b);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"t") == 0) {
    param.setValue(t);
    return param.addObject(3, this);
  }

  return -1;
}

int
HSSSectionIntegration::updateParameter(int parameterID,
					      Information &info)
{
  switch (parameterID) {
  case 1:
    h = info.theDouble;
    return 0;
  case 2:
    b = info.theDouble;
    return 0;
  case 3:
    t = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
HSSSectionIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
HSSSectionIntegration::getLocationsDeriv(int nFibers, double *dyidh, double *dzidh)
{
  for (int i = 0; i < nFibers; i++)
    dyidh[i] = 0.0;
  if (dzidh != 0) {
    for (int i = 0; i < nFibers; i++)
      dzidh[i] = 0.0;
  }

  if (parameterID < 1 || parameterID > 3)
    return;

  double dhdh = 0.0;
  double dbdh = 0.0;
  double dtdh = 0.0;
  if (parameterID == 1)
    dhdh = 1.0;
  if (parameterID == 2)
    dbdh = 1.0;
  if (parameterID == 3)
    dtdh = 1.0;
  
  //double hw = h - 2*t;
  //double bw = b - 2*t;  
  double dhwdh = dhdh - 2*dtdh;
  double dbwdh = dbdh - 2*dtdh;
  
  int i, loc;
  
  //double yIncr  = hw/Nfh;
  //double yStart = 0.5 * (hw-yIncr);
  double dyIncrdh = dhwdh/Nfh;
  double dyStartdh = 0.5*(dhwdh-dyIncrdh);
  
  for (loc = 0, i = 0; loc < Nfh; loc++, i++) {
    //yi[loc]     = yStart - yIncr*i;
    //yi[loc+Nfh] = yi[loc];
    dyidh[loc]     = dyStartdh - dyIncrdh*i;
    dyidh[loc+Nfh] = dyidh[loc];    
  }
  if (dzidh != 0) {
    for (loc = 0; loc < Nfh; loc++) {
      //zi[loc]     = 0.5*(bw + t);
      //zi[loc+Nfh] = -zi[loc];
      dzidh[loc]     = 0.5*(dbwdh + dtdh);
      dzidh[loc+Nfh] = -dzidh[loc];      
    }
  }
  
  for (loc = 2*Nfh, i = 0; i < Nfb; loc++, i++) {
    //yi[loc]     = 0.5*(hw + t);
    //yi[loc+Nfb] = -yi[loc];
    dyidh[loc]     = 0.5*(dhwdh + dtdh);
    dyidh[loc+Nfb] = -dyidh[loc];    
  }
  if (dzidh != 0) {
    //double zIncr  = bw/Nfb;
    //double zStart = 0.5 * (bw-zIncr);
    double dzIncrdh  = dbwdh/Nfb;
    double dzStartdh = 0.5 * (dbwdh-dzIncrdh);    
    for (loc = 2*Nfh, i = 0; i < Nfb; loc++, i++) {
      //zi[loc]     = zStart - zIncr*i;
      //zi[loc+Nfb] = zi[loc];
      dzidh[loc]     = dzStartdh - dzIncrdh*i;
      dzidh[loc+Nfb] = dzidh[loc];      
    }
  }

  // Four corners
  loc = 2*Nfh + 2*Nfb;
  //double yy = 0.5*(hw+t);
  double dyydh = 0.5*(dhwdh+dtdh);  
  //yi[loc++] =  yy;
  //yi[loc++] = -yy;
  //yi[loc++] = -yy;
  //yi[loc++] =  yy;
  dyidh[loc++] =  dyydh;
  dyidh[loc++] = -dyydh;
  dyidh[loc++] = -dyydh;
  dyidh[loc++] =  dyydh;  
  if (dzidh != 0) {
    loc = 2*Nfh + 2*Nfb;
    //double zz = 0.5*(bw+t);
    double dzzdh = 0.5*(dbwdh+dtdh);    
    //zi[loc++] =  zz;
    //zi[loc++] =  zz;
    //zi[loc++] = -zz;
    //zi[loc++] = -zz;
    dzidh[loc++] =  dzzdh;
    dzidh[loc++] =  dzzdh;
    dzidh[loc++] = -dzzdh;
    dzidh[loc++] = -dzzdh;    
  }

  return;
}

void
HSSSectionIntegration::getWeightsDeriv(int nFibers, double *dwtdh)
{
  for (int i = 0; i < nFibers; i++)
    dwtdh[i] = 0.0;

  if (parameterID < 1 || parameterID > 3)
    return;

  double dhdh = 0.0;
  double dbdh = 0.0;
  double dtdh = 0.0;
  if (parameterID == 1)
    dhdh = 1.0;
  if (parameterID == 2)
    dbdh = 1.0;
  if (parameterID == 3)
    dtdh = 1.0;

  double hw = h - 2*t;
  double bw = b - 2*t;
  double dhwdh = dhdh - 2*dtdh;
  double dbwdh = dbdh - 2*dtdh;
  
  // Assuming Nft = 1
  //double a_h = hw*t/(Nfh);
  //double a_b = bw*t/(Nfb);
  double da_hdh = (hw*dtdh + dhwdh*t)/Nfh;
  double da_bdh = (bw*dtdh + dbwdh*t)/Nfb;
  
  int i, loc;
  
  for (loc = 0; loc < Nfh; loc++) {
    //wt[loc]     = a_h;
    //wt[loc+Nfh] = a_h;
    dwtdh[loc]     = da_hdh;
    dwtdh[loc+Nfh] = da_hdh;    
  }

  for (loc = 2*Nfh, i = 0; i < Nfb; i++, loc++) {
    //wt[loc]     = a_b;
    //wt[loc+Nfb] = a_b;
    dwtdh[loc]     = da_bdh;
    dwtdh[loc+Nfb] = da_bdh;        
  }

  // Four corners
  //double a_t = t*t;
  double da_tdh = 2*t*dtdh;
  loc = 2*Nfh + 2*Nfb;
  for (i = 0; i < 4; i++, loc++) {
    //wt[loc] = a_t;
    dwtdh[loc] = da_tdh;
  }
  
  return;
}

void
HSSSectionIntegration::Print(OPS_Stream &s, int flag)
{
  s << "HSS" << endln;
  s << " h = "  << h << endln;
  s << " b = " << b << endln;
  s << " t = " << t << endln;
  s << " Nfh = " << Nfh << endln;
  s << " Nfb = " << Nfb << endln;  
  s << " Nft = " << Nft << endln;

  return;
}

int
HSSSectionIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(6);

  data(0) = h;
  data(1) = b;
  data(2) = t;
  data(3) = Nfh;
  data(4) = Nfb;
  data(5) = Nft;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "HSSSectionIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
HSSSectionIntegration::recvSelf(int cTag, Channel &theChannel,
				FEM_ObjectBroker &theBroker)
{
  static Vector data(6);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "HSSSectionIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  h   = data(0);
  b   = data(1);
  t   = data(2);
  Nfh = (int)data(3);
  Nfb = (int)data(4);
  Nft = (int)data(5);

  return 0;
}
