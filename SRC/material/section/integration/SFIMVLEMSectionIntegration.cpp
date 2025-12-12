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

// $Revision: 1.6 $
// $Date: 2010-09-13 21:31:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/integration/SFIMVLEMSectionIntegration.cpp,v $

#include <SFIMVLEMSectionIntegration.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

#include <elementAPI.h>
#include <NDFiberSection2d.h>
#include <NDMaterial.h>

void* OPS_SFIMVLEMSection2d()
{
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: section SFIMVLEM tag? N? '-thick' *bList '-width' *hList '-mat' *conc" << endln;
    return 0;
  }
  
  // int
  int numdata = 2;
  int idata[2];
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid section SFIMVLEM2d int inputs" << endln;
    return 0;
  }
  
  int tag = idata[0];
  int Nlayers = idata[1];

  int *concTags = new int[Nlayers];
  double *thick = new double[Nlayers];
  double *width = new double[Nlayers];
  
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *str = OPS_GetString();
    int error = 0;
    if (strncmp(str,"-mat",80) == 0) {
      numdata = Nlayers;
      if (OPS_GetIntInput(&numdata, concTags) != 0) {
	opserr << "SFIMVLEMSection error reading concrete tags" << endln;
	error = -1;
      }
    }
    if (strncmp(str,"-thick",80) == 0) {
      numdata = Nlayers;
      if (OPS_GetDoubleInput(&numdata, thick) != 0) {
	opserr << "SFIMVLEMSection error reading thicknesses" << endln;
	error = -1;
      }
    }
    if (strncmp(str,"-width",80) == 0) {
      numdata = Nlayers;
      if (OPS_GetDoubleInput(&numdata, width) != 0) {
	opserr << "SFIMVLEMSection error reading widths" << endln;
	error = -1;
      }
    }

    if (error < 0) {
      delete [] concTags;
      delete [] thick;
      delete [] width;
      
      return 0;      
    }
  }

  NDMaterial **theMats = new NDMaterial *[Nlayers];
  for (int i = 0; i < Nlayers; i++) {
    NDMaterial *concrete = OPS_getNDMaterial(concTags[i]);
    int error = 0;
    if (concrete == 0) {
      opserr << "WARNING ND material does not exist\n";
      opserr << "material: " << concTags[i];
      opserr << "\nSFIMVLEM2d section: " << tag << endln;
      error = -1;
    }
    theMats[i] = concrete;
    
    if (error < 0) {
      delete [] theMats;
      return 0;
    }
  }
  
  SFIMVLEMSectionIntegration rcsect(Nlayers, width, thick);
  
  // Parsing was successful, allocate the section
  SectionForceDeformation* theSection = new NDFiberSection2d(tag, Nlayers, theMats, rcsect);
  
  delete [] theMats;
  return theSection;
}

SFIMVLEMSectionIntegration::SFIMVLEMSectionIntegration(int N,
						       double *H,
						       double *B):
  SectionIntegration(SECTION_INTEGRATION_TAG_RC),
  Nlayers(N), h(0), b(0), parameterID(0)
{
  if (Nlayers > 0) {
    h = new double [Nlayers];
    for (int i = 0; H != 0 && i < Nlayers; i++)
      h[i] = H[i];
    b = new double [Nlayers];    
    for (int i = 0; B != 0 && i < Nlayers; i++)
      b[i] = B[i];
  }
}

SFIMVLEMSectionIntegration::SFIMVLEMSectionIntegration():
  SectionIntegration(SECTION_INTEGRATION_TAG_RC),
  Nlayers(0), h(0), b(0), parameterID(0)
{
  
}

SFIMVLEMSectionIntegration::~SFIMVLEMSectionIntegration()
{
  if (h != 0)
    delete [] h;
  if (b != 0)
    delete [] b;
}

int
SFIMVLEMSectionIntegration::getNumFibers(FiberType type)
{
  if (type == concrete)
    return Nlayers;
  if (type == all)
    return Nlayers;

  return 0;
}

int
SFIMVLEMSectionIntegration::arrangeFibers(NDMaterial **theMaterials,
					  NDMaterial **theReinforcedConcrete)
{
  for (int i = 0; i < Nlayers; i++)
    theMaterials[i] = theReinforcedConcrete[i];

  return 0;
}

void
SFIMVLEMSectionIntegration::getFiberLocations(int nFibers, double *yi, double *zi)
{
  double Lwall = 0.0;
  for (int i = 0; i < Nlayers; i++)
    Lwall += h[i];

  for (int i = 0; i < Nlayers; i++) {
    double sumh = 0.0;
    for (int j = 0; j < i+1; j++)
      sumh += h[j];
    
    yi[i] = (sumh - 0.5*h[i]) - 0.5*Lwall;
    yi[i+Nlayers] = yi[i];
  }

  if (zi != 0) {
    for (int i = 0; i < 2*Nlayers; i++)
      zi[i] = 0.0;
  }

  return;
}

void
SFIMVLEMSectionIntegration::getFiberWeights(int nFibers, double *wt)
{
  for (int i = 0; i < Nlayers; i++)
    wt[i] = b[i]*h[i];

  return;
}

SectionIntegration*
SFIMVLEMSectionIntegration::getCopy(void)
{
  return new SFIMVLEMSectionIntegration(Nlayers, h, b);
}

int
SFIMVLEMSectionIntegration::setParameter(const char **argv, int argc,
				Parameter &param)
{
  if (argc < 1)
    return -1;

  return -1;
}

int
SFIMVLEMSectionIntegration::updateParameter(int parameterID,
				   Information &info)
{
  return -1;
}

int
SFIMVLEMSectionIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
SFIMVLEMSectionIntegration::getLocationsDeriv(int nFibers, double *dyidh, double *dzidh)
{
  return;
}

void
SFIMVLEMSectionIntegration::getWeightsDeriv(int nFibers, double *dwtsdh)
{
  return;
}

void
SFIMVLEMSectionIntegration::Print(OPS_Stream &s, int flag)
{
  s << "SFIMVLEM" << endln;
  s << " Nlayers = " << Nlayers << endln;
  for (int i = 0; i < Nlayers; i++) {
    s << "   " << h[i] << ' ' << b[i] << endln;
  }

  return;
}

int
SFIMVLEMSectionIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(9);

  /*
  data(0) = d;
  data(1) = b;
  data(2) = Atop;
  data(8) = Abottom;
  data(3) = Aside;
  data(4) = cover;
  data(5) = Nfcore;
  data(6) = Nfcover;
  data(7) = Nfs;
  */
  
  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "SFIMVLEMSectionIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
SFIMVLEMSectionIntegration::recvSelf(int cTag, Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
{
  static Vector data(9);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "SFIMVLEMSectionIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }

  /*
  d       = data(0);
  b       = data(1);
  Atop    = data(2);
  Abottom = data(8);
  Aside   = data(3);
  cover   = data(4);
  Nfcore  = (int)data(5);
  Nfcover = (int)data(6);
  Nfs     = (int)data(7);
  */
  
  return 0;
}
