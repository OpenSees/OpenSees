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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/integration/RCWallSectionIntegration.cpp,v $

#include <RCWallSectionIntegration.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

#include <elementAPI.h>
#include <FiberSection2d.h>
#include <UniaxialMaterial.h>

void* OPS_RCWall2d()
{
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING insufficient arguments\n";
    // '-thick',*bList,'-width',*hList,'-rho',*rhoList,'-matConcrete',*conc,'-matSteel',*steel,
    opserr << "Want: section RCWall2d tag? N? *concTag? *steelTag? *h *b *rho" << endln;
    return 0;
  }
  
  // int
  int numdata = 2;
  int idata[2];
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid section RCWall2d int inputs" << endln;
    return 0;
  }
  
  int tag = idata[0];
  int Nlayers = idata[1];

  int *concTags = new int[Nlayers];
  int *steelTags = new int[Nlayers];
  double *thick = new double[Nlayers];
  double *width = new double[Nlayers];
  double *rho = new double[Nlayers];  
  
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *str = OPS_GetString();
    int error = 0;
    if (strncmp(str,"-matConcrete",80) == 0) {
      numdata = Nlayers;
      if (OPS_GetIntInput(&numdata, concTags) != 0) {
	opserr << "RCWallSection error reading concrete tags" << endln;
	error = -1;
      }
    }
    if (strncmp(str,"-matSteel",80) == 0) {
      numdata = Nlayers;
      if (OPS_GetIntInput(&numdata, steelTags) != 0) {
	opserr << "RCWallSection error reading steel tags" << endln;
	error = -1;
      }
    }
    if (strncmp(str,"-thick",80) == 0) {
      numdata = Nlayers;
      if (OPS_GetDoubleInput(&numdata, thick) != 0) {
	opserr << "RCWallSection error reading thicknesses" << endln;
	error = -1;
      }
    }
    if (strncmp(str,"-width",80) == 0) {
      numdata = Nlayers;
      if (OPS_GetDoubleInput(&numdata, width) != 0) {
	opserr << "RCWallSection error reading widths" << endln;
	error = -1;
      }
    }
    if (strncmp(str,"-rho",80) == 0) {
      numdata = Nlayers;
      if (OPS_GetDoubleInput(&numdata, rho) != 0) {
	opserr << "RCWallSection error reading rho values" << endln;
	error = -1;
      }
    }        

    if (error < 0) {
      delete [] concTags;
      delete [] steelTags;
      delete [] thick;
      delete [] width;
      delete [] rho;
      
      return 0;      
    }
  }

  UniaxialMaterial **theMats = new UniaxialMaterial *[2*Nlayers];
  for (int i = 0; i < Nlayers; i++) {
    UniaxialMaterial *concrete = OPS_getUniaxialMaterial(concTags[i]);
    int error = 0;
    if (concrete == 0) {
      opserr << "WARNING uniaxial material does not exist\n";
      opserr << "material: " << concTags[i];
      opserr << "\nRCWall2d section: " << tag << endln;
      error = -1;
    }
    theMats[i] = concrete;
    
    UniaxialMaterial *steel = OPS_getUniaxialMaterial(steelTags[i]);
    if (steel == 0) {
      opserr << "WARNING uniaxial material does not exist\n";
      opserr << "material: " << steelTags[i];
      opserr << "\nRCWall2d section: " << tag << endln;
      error = -1;
    }
    theMats[i+Nlayers] = steel;
    
    if (error < 0) {
      delete [] theMats;
      return 0;
    }
  }
  
  RCWallSectionIntegration rcsect(Nlayers, width, thick, rho);
  
  // Parsing was successful, allocate the section
  SectionForceDeformation* theSection = new FiberSection2d(tag, 2*Nlayers, theMats, rcsect);
  
  delete [] theMats;
  return theSection;
}

RCWallSectionIntegration::RCWallSectionIntegration(int N,
				     double *H,
				     double *B,
				     double *R):
  SectionIntegration(SECTION_INTEGRATION_TAG_RC),
  Nlayers(N), h(0), b(0), rho(0), parameterID(0)
{
  if (Nlayers > 0) {
    h = new double [Nlayers];
    for (int i = 0; H != 0 && i < Nlayers; i++)
      h[i] = H[i];
    b = new double [Nlayers];    
    for (int i = 0; B != 0 && i < Nlayers; i++)
      b[i] = B[i];
    rho = new double [Nlayers];    
    for (int i = 0; R != 0 && i < Nlayers; i++)
      rho[i] = R[i];    
  }
}

RCWallSectionIntegration::RCWallSectionIntegration():
  SectionIntegration(SECTION_INTEGRATION_TAG_RC),
  Nlayers(0), h(0), b(0), rho(0), parameterID(0)
{
  
}

RCWallSectionIntegration::~RCWallSectionIntegration()
{
  if (h != 0)
    delete [] h;
  if (b != 0)
    delete [] b;
  if (rho != 0)
    delete [] rho;  
}

int
RCWallSectionIntegration::getNumFibers(FiberType type)
{
  if (type == steel)
    return Nlayers;
  if (type == concrete)
    return Nlayers;
  if (type == all)
    return 2*Nlayers;

  return 0;
}

int
RCWallSectionIntegration::arrangeFibers(UniaxialMaterial **theMaterials,
				    UniaxialMaterial *theConcrete,
				    UniaxialMaterial *theSteel)
{
  int i;
  for (i = 0; i < Nlayers; i++)
    theMaterials[i] = theConcrete;
  for ( ; i < 2*Nlayers; i++)
    theMaterials[i] = theSteel;

  return 0;
}

void
RCWallSectionIntegration::getFiberLocations(int nFibers, double *yi, double *zi)
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
RCWallSectionIntegration::getFiberWeights(int nFibers, double *wt)
{
  for (int i = 0; i < Nlayers; i++) {
    double As = rho[i]*b[i]*h[i];
    wt[i+Nlayers] = As;
    wt[i] = b[i]*h[i] - As;
  }

  return;
}

SectionIntegration*
RCWallSectionIntegration::getCopy(void)
{
  return new RCWallSectionIntegration(Nlayers, h, b, rho);
}

int
RCWallSectionIntegration::setParameter(const char **argv, int argc,
				Parameter &param)
{
  if (argc < 1)
    return -1;

  return -1;
}

int
RCWallSectionIntegration::updateParameter(int parameterID,
				   Information &info)
{
  return -1;
}

int
RCWallSectionIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
RCWallSectionIntegration::getLocationsDeriv(int nFibers, double *dyidh, double *dzidh)
{
  return;
}

void
RCWallSectionIntegration::getWeightsDeriv(int nFibers, double *dwtsdh)
{
  return;
}

void
RCWallSectionIntegration::Print(OPS_Stream &s, int flag)
{
  s << "RC Wall" << endln;
  s << " Nlayers = " << Nlayers << endln;
  for (int i = 0; i < Nlayers; i++) {
    s << "   " << h[i] << ' ' << b[i] << ' ' << rho[i] << endln;
  }

  return;
}

int
RCWallSectionIntegration::sendSelf(int cTag, Channel &theChannel)
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
    opserr << "RCWallSectionIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
RCWallSectionIntegration::recvSelf(int cTag, Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
{
  static Vector data(9);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "RCWallSectionIntegration::recvSelf() - failed to receive Vector data\n";
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
