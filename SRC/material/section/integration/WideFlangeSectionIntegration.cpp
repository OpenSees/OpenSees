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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/integration/WideFlangeSectionIntegration.cpp,v $

#include <WideFlangeSectionIntegration.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

#include <elementAPI.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <FiberSection2d.h>
#include <NDFiberSection2d.h>
#include <NDFiberSectionWarping2d.h>

void * OPS_ADD_RUNTIME_VPV(OPS_WFSection2d)
{
  if (OPS_GetNumRemainingInputArgs() < 8) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: section WFSection2d tag? matTag? d? tw? bf? tf? nfdw? nftf? <-nd shape?>" << endln;
    return 0;
  }
  
  int tag, matTag;
  double d, tw, bf, tf;
  int nfdw, nftf;
  
  SectionForceDeformation* theSection = 0;
  
  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid section WFSection2d tag" << endln;
    return 0;
  }
  
  if (OPS_GetIntInput(&numdata, &matTag) < 0) {
    opserr << "WARNING invalid section WFSection2d matTag" << endln;
    return 0;
  }
  
  if (OPS_GetDoubleInput(&numdata, &d) < 0) {
    opserr << "WARNING invalid d" << endln;
    opserr << "WFSection2d section: " << tag << endln;
    return 0;
  }
  
  if (OPS_GetDoubleInput(&numdata, &tw) < 0) {
    opserr << "WARNING invalid tw" << endln;
    opserr << "WFSection2d section: " << tag << endln;
    return 0;
  }
  
  if (OPS_GetDoubleInput(&numdata, &bf) < 0) {
    opserr << "WARNING invalid bf" << endln;
    opserr << "WFSection2d section: " << tag << endln;
    return 0;
  }
  
  if (OPS_GetDoubleInput(&numdata, &tf) < 0) {
    opserr << "WARNING invalid tf" << endln;
    opserr << "WFSection2d section: " << tag << endln;
    return 0;
  }
  
  if (OPS_GetIntInput(&numdata, &nfdw) < 0) {
    opserr << "WARNING invalid nfdw" << endln;
    opserr << "WFSection2d section: " << tag << endln;
    return 0;
  }
  
  if (OPS_GetIntInput(&numdata, &nftf) < 0) {
    opserr << "WARNING invalid nftf" << endln;
    opserr << "WFSection2d section: " << tag << endln;
    return 0;
  }
  
  WideFlangeSectionIntegration wfsect(d, tw, bf, tf, nfdw, nftf);
  
  int numFibers = wfsect.getNumFibers();
  
  if (OPS_GetNumRemainingInputArgs() > 0) {
    
    double shape = 1.0;
    if (OPS_GetNumRemainingInputArgs() > 1) {
      if (OPS_GetDoubleInput(&numdata, &shape) < 0) {
	opserr << "WARNING invalid shape" << endln;
	opserr << "WFSection2d section: " << tag << endln;
	return 0;
      }
    }
    
    NDMaterial *theSteel = OPS_getNDMaterial(matTag);
    
    if (theSteel == 0) {
      opserr << "WARNING ND material does not exist\n";
      opserr << "material: " << matTag;
      opserr << "\nWFSection2d section: " << tag << endln;
      return 0;
    }
    
    NDMaterial **theMats = new NDMaterial *[numFibers];
    
    wfsect.arrangeFibers(theMats, theSteel);
    
    // Parsing was successful, allocate the section
    theSection = 0;
    if (OPS_GetNumRemainingInputArgs() > 0) {
      const char* flag = OPS_GetString();
      if (strcmp(flag,"-nd") == 0) {
	theSection = new NDFiberSection2d(tag, numFibers, theMats, wfsect, shape);
      } else if (strcmp(flag,"-ndWarping") == 0) {
	theSection = new NDFiberSectionWarping2d(tag, numFibers, theMats, wfsect, shape);
      }
    }
    delete [] theMats;
  }
  else {
    UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(matTag);
    
    if (theSteel == 0) {
      opserr << "WARNING uniaxial material does not exist\n";
      opserr << "material: " << matTag;
      opserr << "\nWFSection2d section: " << tag << endln;
      return 0;
    }
    
    UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];
    
    wfsect.arrangeFibers(theMats, theSteel);
    
    // Parsing was successful, allocate the section
    theSection = new FiberSection2d(tag, numFibers, theMats, wfsect);
    
    delete [] theMats;
  }
  
  return theSection;
}  

WideFlangeSectionIntegration::WideFlangeSectionIntegration(double D,
							   double TW,
							   double BF,
							   double TF,
							   int NFDW,
							   int NFTF):
  SectionIntegration(SECTION_INTEGRATION_TAG_WideFlange),
  d(D), tw(TW), bf(BF), tf(TF), Nfdw(NFDW), Nftf(NFTF), parameterID(0)
{
  
}

WideFlangeSectionIntegration::WideFlangeSectionIntegration():
  SectionIntegration(SECTION_INTEGRATION_TAG_WideFlange),
  d(0.0), tw(0.0), bf(0.0), tf(0.0), Nfdw(0), Nftf(0), parameterID(0)
{
  
}

WideFlangeSectionIntegration::~WideFlangeSectionIntegration()
{
  
}

int
WideFlangeSectionIntegration::getNumFibers(FiberType type)
{
  return Nfdw + 2*Nftf;
}

int
WideFlangeSectionIntegration::arrangeFibers(UniaxialMaterial **theMaterials,
					    UniaxialMaterial *theSteel)
{
  int numFibers = this->getNumFibers();

  for (int i = 0; i < numFibers; i++)
    theMaterials[i] = theSteel;
  
  return 0;
}

int
WideFlangeSectionIntegration::arrangeFibers(NDMaterial **theMaterials,
					    NDMaterial *theSteel)
{
  int numFibers = this->getNumFibers();

  for (int i = 0; i < numFibers; i++)
    theMaterials[i] = theSteel;
  
  return 0;
}

void
WideFlangeSectionIntegration::getFiberLocations(int nFibers, double *yi, double *zi)
{
  double dw = d-2*tf;
  
  int loc;
  
  double yIncr  = tf/Nftf;
  double yStart = 0.5 * (d-yIncr);
  
  for (loc = 0; loc < Nftf; loc++) {
    yi[loc] = yStart - yIncr*loc;
    yi[nFibers-loc-1] = -yi[loc];
  }
  
  yIncr  = dw/Nfdw;
  yStart = 0.5 * (dw-yIncr);
  
  for (int count = 0; loc < nFibers-Nftf; loc++, count++) {
    yi[loc] = yStart - yIncr*count;
  }

  if (zi != 0) {
    for (int i = 0; i < nFibers; i++)
      zi[i] = 0.0;
  }

  return;
}

void
WideFlangeSectionIntegration::getFiberWeights(int nFibers, double *wt)
{
  double dw = d-2*tf;
  
  double a_f = bf*tf/Nftf;
  double a_w = dw*tw/Nfdw;
  
  int loc = 0;
  
  for (loc = 0; loc < Nftf; loc++) {
    wt[loc] = a_f;
    wt[nFibers-loc-1] = a_f;
  }
  
  for ( ; loc < nFibers-Nftf; loc++) {
    wt[loc] = a_w;
  }

  return;
}

SectionIntegration*
WideFlangeSectionIntegration::getCopy(void)
{
  return new WideFlangeSectionIntegration(d, tw, bf, tf, Nfdw, Nftf);
}

int
WideFlangeSectionIntegration::setParameter(const char **argv, int argc,
					   Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"d") == 0 || strcmp(argv[0],"db") == 0) {
    param.setValue(d);
    return param.addObject(1, this);    
  }
  if (strcmp(argv[0],"tw") == 0) {
    param.setValue(tw);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"bf") == 0 || strcmp(argv[0],"b") == 0) {
    param.setValue(bf);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"tf") == 0) {
    param.setValue(tf);
    return param.addObject(4, this);
  }

  return -1;
}

int
WideFlangeSectionIntegration::updateParameter(int parameterID,
					      Information &info)
{
  switch (parameterID) {
  case 1:
    d = info.theDouble;
    return 0;
  case 2:
    tw = info.theDouble;
    return 0;
  case 3:
    bf = info.theDouble;
    return 0;
  case 4:
    tf = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
WideFlangeSectionIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
WideFlangeSectionIntegration::getLocationsDeriv(int nFibers, double *dyidh, double *dzidh)
{
  //double dw = d-2*tf;
  
  double dddh  = 0.0;
  double ddwdh = 0.0;
  //double dtwdh = 0.0;
  //double dbfdh = 0.0;
  double dtfdh = 0.0;

  if (parameterID == 1) { // d
    dddh  = 1.0;
    ddwdh = 1.0;
  }
  //if (parameterID == 2) // tw
  //  dtwdh = 1.0;
  //if (parameterID == 3) // bf
  //  dbfdh = 1.0;
  if (parameterID == 4) { // tf
    dtfdh =  1.0;
    ddwdh = -2.0;
  }

  //double yIncr  = tf/Nftf;
  //double yStart = 0.5*d - 0.5*yIncr;
  double dyIncrdh  = dtfdh/Nftf;
  double dyStartdh = 0.5 * (dddh-dyIncrdh);

  int loc;
  
  for (loc = 0; loc < Nftf; loc++) {
    //yi[loc] = yStart - yIncr*loc;
    //yi[nFibers-loc-1] = -yi[loc];
    dyidh[loc] = dyStartdh - dyIncrdh*loc;
    dyidh[nFibers-loc-1] = -dyidh[loc];
  }
  
  //yIncr  = dw/Nfdw;
  //yStart = 0.5*dw - 0.5*yIncr;
  dyIncrdh  = ddwdh/Nfdw;
  dyStartdh = 0.5 * (ddwdh-dyIncrdh);
  
  for (int count = 0; loc < nFibers-Nftf; loc++, count++) {
    //yi[loc] = yStart - yIncr*count;
    dyidh[loc] = dyStartdh - dyIncrdh*count;
  }

  if (dzidh != 0) {
    for (int i = 0; i < nFibers; i++)
      dzidh[i] = 0.0;
  }

  return;
}

void
WideFlangeSectionIntegration::getWeightsDeriv(int nFibers, double *dwtdh)
{
  double dw = d-2*tf;
  
  double ddwdh = 0.0;
  double dtwdh = 0.0;
  double dbfdh = 0.0;
  double dtfdh = 0.0;

  if (parameterID == 1) // d
    ddwdh = 1.0;
  if (parameterID == 2) // tw
    dtwdh = 1.0;
  if (parameterID == 3) // bf
    dbfdh = 1.0;
  if (parameterID == 4) { // tf
    dtfdh =  1.0;
    ddwdh = -2.0;
  }
    
  double dAfdh = (bf*dtfdh + dbfdh*tf) / Nftf;
  double dAwdh = (dw*dtwdh + ddwdh*tw) / Nfdw;

  int loc = 0;
  
  for (loc = 0; loc < Nftf; loc++) {
    dwtdh[loc] = dAfdh;
    dwtdh[nFibers-loc-1] = dAfdh;
  }
  
  for ( ; loc < nFibers-Nftf; loc++)
    dwtdh[loc] = dAwdh;

  //for (int i = 0; i < nFibers; i++)
  //  opserr << dwtsdh[i] << ' ';
  //opserr << endln;

  return;
}

void
WideFlangeSectionIntegration::Print(OPS_Stream &s, int flag)
{
  s << "WideFlange" << endln;
  s << " d = "  << d;
  s << " tw = " << tw; 
  s << " bf = " << bf;
  s << " tf = " << tf << endln;
  s << " Nfdw = " << Nfdw;
  s << " Nftf = " << Nftf << endln;

  return;
}

int
WideFlangeSectionIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(6);

  data(0) = d;
  data(1) = tw;
  data(2) = bf;
  data(3) = tf;
  data(4) = Nfdw;
  data(5) = Nftf;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "WideFlangeSectionIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
WideFlangeSectionIntegration::recvSelf(int cTag, Channel &theChannel,
				       FEM_ObjectBroker &theBroker)
{
  static Vector data(6);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "WideFlangeSectionIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  d    = data(0);
  tw   = data(1);
  bf   = data(2);
  tf   = data(3);
  Nfdw = (int)data(4);
  Nftf = (int)data(5);

  return 0;
}
