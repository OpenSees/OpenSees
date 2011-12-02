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

// $Revision: 1.1 $
// $Date: 2006-08-11 18:32:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/integration/RCSectionIntegration.cpp,v $

#include <RCSectionIntegration.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

RCSectionIntegration::RCSectionIntegration(double D,
					   double B,
					   double AM,
					   double AS,
					   double COV,
					   int NFCORE,
					   int NFCOVER,
					   int NFS):
  SectionIntegration(SECTION_INTEGRATION_TAG_RC),
  d(D), b(B), Amain(AM), Aside(AS), cover(COV),
  Nfcore(NFCORE), Nfcover(NFCOVER), Nfs(NFS), parameterID(0)
{
  if (Nfcore < 1)
    Nfcore = 1;

  if (Nfcover < 1)
    Nfcover = 1;

  if (Nfs < 2)
    Nfs = 2;
}

RCSectionIntegration::RCSectionIntegration():
  SectionIntegration(SECTION_INTEGRATION_TAG_RC),
  d(0.0), b(0.0), Amain(0.0), Aside(0.0), cover(0.0),
  Nfcore(1), Nfcover(1), Nfs(2), parameterID(0)
{
  
}

RCSectionIntegration::~RCSectionIntegration()
{
  
}

int
RCSectionIntegration::getNumFibers(void)
{
  return 2*(Nfcore+Nfcover) + Nfs;
}

int
RCSectionIntegration::arrangeFibers(UniaxialMaterial **theMaterials,
				    UniaxialMaterial *theCore,
				    UniaxialMaterial *theCover,
				    UniaxialMaterial *theSteel)
{
  int numFibers = this->getNumFibers();

  int i;
  for (i = 0; i < Nfcore; i++)
    theMaterials[i] = theCore;
  for ( ; i < numFibers-Nfs; i++)
    theMaterials[i] = theCover;
  for ( ; i < numFibers; i++)
    theMaterials[i] = theSteel;

  return 0;
}

void
RCSectionIntegration::getFiberLocations(int nFibers, double *xi)
{
  int loc;
  int i;

  double yIncr  = (d-2*cover)/Nfcore;
  double yStart = 0.5 * ((d-2*cover)-yIncr);
  
  for (loc = 0; loc < Nfcore; loc++) {
    xi[loc] = yStart - yIncr*loc;
    xi[loc+Nfcore] = xi[loc];
  }

  loc += Nfcore;
  
  yIncr = cover/Nfcover;
  yStart = 0.5 * (d-yIncr);

  for (i = 0; i < Nfcover; i++, loc++) {
    xi[loc] = yStart - yIncr*i;
    xi[loc+Nfcover] = -xi[loc];
  }

  loc += Nfcover;

  xi[loc++] =  0.5*d-cover;
  xi[loc++] = -0.5*d+cover;

  if (Nfs > 2) {
    double spacing = (d-2*cover)/(Nfs-1);
    for (int i = 1; i <= Nfs-2; i++)
      xi[loc++] = (-0.5*d+cover) + spacing*i;
  }

  return;
}

void
RCSectionIntegration::getFiberWeights(int nFibers, double *wt)
{
  int loc;
  int i;

  double Acore  = (b-2*cover)*(d-2*cover)/Nfcore;
  double Acover = (2*cover)*(d-2*cover)/Nfcore;
    
  for (loc = 0; loc < Nfcore; loc++)
    wt[loc] = Acore;

  for (i = 0; i < Nfcore; i++, loc++)
    wt[loc] = Acover;

  Acover = (cover*b)/Nfcover;

  for (i = 0; i < 2*Nfcover; i++, loc++)
    wt[loc] = Acover;
  
  wt[loc++] = Nfs*Amain;
  wt[loc++] = Nfs*Amain;

  for ( ; loc < nFibers; loc++)
    wt[loc] = 2*Aside;

  return;
}

SectionIntegration*
RCSectionIntegration::getCopy(void)
{
  return new RCSectionIntegration(d, b, Amain, Aside, cover,
				  Nfcore, Nfcover, Nfs);
}

int
RCSectionIntegration::setParameter(const char **argv, int argc,
				   Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"d") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"b") == 0)
    return param.addObject(2, this);

  if (strcmp(argv[0],"Amain") == 0)
    return param.addObject(3, this);

  if (strcmp(argv[0],"Aside") == 0)
    return param.addObject(4, this);

  if (strcmp(argv[0],"As") == 0)
    return param.addObject(5, this);

  if (strcmp(argv[0],"cover") == 0)
    return param.addObject(6, this);

  return -1;
}

int
RCSectionIntegration::updateParameter(int parameterID,
				      Information &info)
{
  switch (parameterID) {
  case 1:
    d = info.theDouble;
    return 0;
  case 2:
    b = info.theDouble;
    return 0;
  case 3:
    Amain = info.theDouble;
    return 0;
  case 4:
    Aside = info.theDouble;
    return 0;
  case 5:
    Amain = Aside = info.theDouble;
    return 0;
  case 6:
    cover = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
RCSectionIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
RCSectionIntegration::getLocationsDeriv(int nFibers, double *dptsdh)
{
  double dddh  = 0.0;
  double dcoverdh = 0.0;
  
  if (parameterID == 1) // d
    dddh  = 1.0;
  if (parameterID == 6) // cover
    dcoverdh =  1.0;

  int loc;
  int i;

  double dyIncrdh  = (dddh-2*dcoverdh)/Nfcore;
  double dyStartdh = 0.5 * ((dddh-2*dcoverdh)-dyIncrdh);
  
  for (loc = 0; loc < Nfcore; loc++) {
    dptsdh[loc] = dyStartdh - dyIncrdh*loc;
    dptsdh[loc+Nfcore] = dptsdh[loc];
  }

  loc += Nfcore;
  
  dyIncrdh = dcoverdh/Nfcover;
  dyStartdh = 0.5 * (dddh-dyIncrdh);

  for (i = 0; i < Nfcover; i++, loc++) {
    dptsdh[loc] = dyStartdh - dyIncrdh*i;
    dptsdh[loc+Nfcover] = -dptsdh[loc];
  }

  loc += Nfcover;

  dptsdh[loc++] =  0.5*dddh-dcoverdh;
  dptsdh[loc++] = -0.5*dddh+dcoverdh;

  if (Nfs > 2) {
    double dspacingdh = (dddh-2*dcoverdh)/(Nfs-1);
    for (int i = 1; i <= Nfs-2; i++)
      dptsdh[loc++] = (-0.5*dddh+dcoverdh) + dspacingdh*i;
  }

  //for (int i = 0; i < nFibers; i++)
  //  opserr << dptsdh[i] << ' ';
  //opserr << endln;

  return;
}

void
RCSectionIntegration::getWeightsDeriv(int nFibers, double *dwtsdh)
{
  double dddh  = 0.0;
  double dbdh = 0.0;
  double dAmaindh = 0.0;
  double dAsidedh = 0.0;
  double dcoverdh = 0.0;
  
  if (parameterID == 1) // d
    dddh  = 1.0;
  if (parameterID == 2) // b
    dbdh = 1.0;
  if (parameterID == 3) // Amain
    dAmaindh = 1.0;
  if (parameterID == 4) // Aside
    dAsidedh = 1.0;
  if (parameterID == 5) // Amain and Aside
    dAmaindh = dAsidedh = 1.0;
  if (parameterID == 6) // cover
    dcoverdh =  1.0;
    
  int loc;
  int i;

  double dAcoredh = ((b-2*cover)*(dddh-2*dcoverdh) +
		     (dbdh-2*dcoverdh)*(d-2*cover)) / Nfcore;
  double dAcoverdh = ((2*cover)*(dddh-2*dcoverdh) + 
		      (2*dcoverdh)*(d-2*cover)) / Nfcore;
    
  for (loc = 0; loc < Nfcore; loc++)
    dwtsdh[loc] = dAcoredh;

  for (i = 0; i < Nfcore; i++, loc++)
    dwtsdh[loc] = dAcoverdh;

  dAcoverdh = (cover*dbdh+dcoverdh*b)/Nfcover;

  for (i = 0; i < 2*Nfcover; i++, loc++)
    dwtsdh[loc] = dAcoverdh;
  
  dwtsdh[loc++] = Nfs*dAmaindh;
  dwtsdh[loc++] = Nfs*dAmaindh;

  for ( ; loc < nFibers; loc++)
    dwtsdh[loc] = 2*dAsidedh;

  return;

  //for (int i = 0; i < nFibers; i++)
  //  opserr << dwtsdh[i] << ' ';
  //opserr << endln;

  return;
}

void
RCSectionIntegration::Print(OPS_Stream &s, int flag)
{
  s << "RC" << endln;
  s << " d = "  << d;
  s << " b = " << b; 
  s << " Amain = " << Amain;
  s << " Amain = " << Aside;
  s << " cover = " << cover << endln;
  s << " Nfcore = " << Nfcore;
  s << " Nfcover = " << Nfcover;
  s << " Nfs = " << Nfs << endln;

  return;
}

int
RCSectionIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(8);

  data(0) = d;
  data(1) = b;
  data(2) = Amain;
  data(3) = Aside;
  data(4) = cover;
  data(5) = Nfcore;
  data(6) = Nfcover;
  data(7) = Nfs;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "RCSectionIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
RCSectionIntegration::recvSelf(int cTag, Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
{
  static Vector data(8);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "RCSectionIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  d       = data(0);
  b       = data(1);
  Amain   = data(2);
  Aside   = data(3);
  cover   = data(4);
  Nfcore  = (int)data(5);
  Nfcover = (int)data(6);
  Nfs     = (int)data(7);

  return 0;
}
