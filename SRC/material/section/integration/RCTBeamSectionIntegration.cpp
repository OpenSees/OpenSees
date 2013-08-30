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

// $Revision: 1.5 $
// $Date: 2007/12/01 01:03:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/integration/RCSectionIntegration.cpp,v $

#include <RCTBeamSectionIntegration.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>

RCTBeamSectionIntegration::RCTBeamSectionIntegration(double D,
						     double BW,
						     double BEFF,
						     double HF,
						     double AT,
						     double AB,
						     double FLCOV,
						     double WCOV,
						     int NFLCOVER,
						     int NWCOVER,
						     int NFLCORE,
						     int NWCORE,
						     int NSTEELTOP,
						     int NSTEELBOTTOM):
  SectionIntegration(SECTION_INTEGRATION_TAG_RCT),
  d(D), bw(BW), beff(BEFF), hf(HF), Atop(AT), Abottom(AB), flcov(FLCOV),
  wcov(WCOV), Nflcover(NFLCOVER), Nwcover(NWCOVER),
  Nflcore(NFLCORE), Nwcore(NWCORE), NsteelTop(NSTEELTOP) , NsteelBottom(NSTEELBOTTOM), parameterID(0)
{
  /*
  if (Nflcover < 1)
    Nflcover = 1;

  if (Nwcover < 1)
    Nwcover = 1;

  if (Nflcore < 1)
    Nflcore = 1;

  if (Nwcore < 1)
    Nwcore = 1;

  if (NsteelTop < 1)
    NsteelTop = 1;
	
  if (NsteelBottom < 1)
    NsteelBottom = 1;
  */
}

RCTBeamSectionIntegration::RCTBeamSectionIntegration():
  SectionIntegration(SECTION_INTEGRATION_TAG_RCT),
  d(0.0), bw(0.0), beff(0.0), hf(0.0), Atop(0.0), Abottom(0.0), flcov(0.0),
  wcov(0.0), Nflcover(0), Nwcover(0), Nflcore(0), Nwcore(0), NsteelTop(0), NsteelBottom(0), parameterID(0)
{
  
}

RCTBeamSectionIntegration::~RCTBeamSectionIntegration()
{
  
}

int
RCTBeamSectionIntegration::getNumFibers(FiberType type)
{
  if (type == steel)
    return NsteelTop + NsteelBottom;
  if (type == concrete)
    return (Nflcover+Nwcover) + (Nflcore+Nwcore);
  if (type == all)
    return (Nflcover+Nwcover) + (Nflcore+Nwcore) + NsteelTop + NsteelBottom;

  return 0;
}

//Check for concreteMcft
int
RCTBeamSectionIntegration::arrangeFibers(UniaxialMaterial **theUni,
					 NDMaterial **theND,
					 NDMaterial *theCore,
					 NDMaterial *theCover,
					 UniaxialMaterial *theSteel)
{
  int loc = 0;

  int Ncore = Nflcore+Nwcore;
  if (theCore != 0)
    for (int i = 0; i < Ncore; i++)
      theND[loc++] = theCore;
  else
    loc += Ncore;

  int Ncover = Nflcover+Nwcover;
  if (theCover != 0) 
    for (int i = 0; i < Ncover; i++)
      theND[loc++] = theCover;

  int Nsteel = NsteelTop+NsteelBottom;
  if (theSteel != 0)
    for (int i = 0; i < Nsteel; i++)
      theUni[i] = theSteel;

  return 0;
}

int
RCTBeamSectionIntegration::arrangeFibers(UniaxialMaterial **theUni,
					 UniaxialMaterial *theCore,
					 UniaxialMaterial *theCover,
					 UniaxialMaterial *theSteel)
{
  int loc = 0;

  int Ncore = Nflcore+Nwcore;
  if (theCore != 0)
    for (int i = 0; i < Ncore; i++)
      theUni[loc++] = theCore;
  else
    loc += Ncore;

  int Ncover = Nflcover+Nwcover;
  if (theCover != 0) 
    for (int i = 0; i < Ncover; i++)
      theUni[loc++] = theCover;
  else 
    loc += Ncover;
  
  int Nsteel = NsteelTop+NsteelBottom;
  if (theSteel != 0)
    for (int i = 0; i < Nsteel; i++)
      theUni[loc++] = theSteel;

  return 0;
}

void
RCTBeamSectionIntegration::getFiberLocations(int nFibers, double *yi, double *zi)
{
  double Yg = ((d-hf)*bw*(d-hf)/2+beff*hf*(d-hf/2)) / (beff*hf+(d-hf)*bw);
  Yg = 0;
  int loc = 0;

  // Locate core web fibers
  double yIncrw  = (d-hf-wcov)/Nwcore;
  double yStartw = -Yg + wcov + 0.5 * yIncrw;  
  for (int i = 0; i < Nwcore; i++)
    yi[loc++] = yStartw + yIncrw*i;
   
  // Locate core flange fibers
  double yIncrfl = (hf-flcov) / Nflcore;
  double yStartfl = d - Yg - hf + yIncrfl*0.5;
  for (int i = 0; i < Nflcore; i++)
    yi[loc++] = yStartfl + yIncrfl*i;
  
  // Locate cover web fibers
  double yIncrwcov = wcov / Nwcover;
  double yStartwcov = -Yg + 0.5 *yIncrwcov;
  for (int i = 0;  i < Nwcover; i++)
    yi[loc++] = yStartwcov + yIncrwcov*i;

  // Locate cover flange fibers
  double yIncrflcov = flcov /Nflcover;
  double yStartflcov = d - Yg - flcov + 0.5*yIncrflcov;
  for (int i = 0; i < Nflcover; i++)
    yi[loc++] = yStartflcov + yIncrflcov*i;

  // Locate Steel fibers
  double yTop = d - Yg - flcov;
  for (int i = 0; i < NsteelTop; i++)
    yi[loc++] = yTop;

  double yBot = -Yg + wcov;
  for (int i = 0; i < NsteelBottom; i++)
    yi[loc++] = yBot;

  if (zi != 0) {
    for (int i = 0; i < nFibers; i++)
      zi[i] = 0.0;
  }
 
  return;
}

void
RCTBeamSectionIntegration::getFiberWeights(int nFibers, double *wt)
{
  int loc = 0;

  // Weight core web fibers
  double Awcore = bw*(d-hf -wcov)/Nwcore;
  for (int i = 0; i < Nwcore; i++)
    wt[loc++] = Awcore;

  // Weight core flange fibers
  double Aflcore = (hf-flcov)*beff / Nflcore;
  for (int i = 0; i < Nflcore; i++)
    wt[loc++] = Aflcore;

  // Weight cover web fibers
  double Awcover = bw * wcov / Nwcover;
  for (int i = 0; i < Nwcover; i++)
    wt[loc++] = Awcover;

  // Weight cover flange fibers
  double Aflcover = beff * flcov/Nflcover;
  for (int i = 0; i < Nflcover; i++)
    wt[loc++] = Aflcover;

  // Weight Steel fibers
  for (int i = 0; i < NsteelTop; i++)
    wt[loc++] = Atop;

  for (int i = 0; i < NsteelBottom; i++)
    wt[loc++] = Abottom;

  return;
}

SectionIntegration*
RCTBeamSectionIntegration::getCopy(void)
{
  return new RCTBeamSectionIntegration(d, bw, beff, hf, Atop, Abottom, 
				       flcov, wcov, Nflcover,
				       Nwcover, Nflcore, Nwcore, 
				       NsteelTop, NsteelBottom);
}

int
RCTBeamSectionIntegration::setParameter(const char **argv, int argc,
				   Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"d") == 0) {
    param.setValue(d);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"bw") == 0) {
    param.setValue(bw);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"beff") == 0) {
    param.setValue(beff);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"hf") == 0) {
    param.setValue(hf);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"Atop") == 0) {
    param.setValue(Atop);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"Abottom") == 0) {
    param.setValue(Abottom);
    return param.addObject(6, this);
  }
  if (strcmp(argv[0],"flcov") == 0) {
    param.setValue(flcov);
    return param.addObject(7, this);
  }
  if (strcmp(argv[0],"wcov") == 0) {
    param.setValue(wcov);
    return param.addObject(8, this);
  }
  return -1;
}

int
RCTBeamSectionIntegration::updateParameter(int parameterID,
				      Information &info)
{
  switch (parameterID) {
  case 1:
    d = info.theDouble;
    return 0;
  case 2:
    bw = info.theDouble;
    return 0;
  case 3:
    beff = info.theDouble;
    return 0;
  case 4:
    hf = info.theDouble;
    return 0;
  case 5:
    Atop = info.theDouble;
    return 0;
  case 6:
    Abottom = info.theDouble;
    return 0;
  case 7:
    flcov = info.theDouble;
    return 0;
  case 8:
    wcov = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
RCTBeamSectionIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
RCTBeamSectionIntegration::getLocationsDeriv(int nFibers, 
					     double *dyidh, 
					     double *dzidh)
{
  double dddh  = 0.0;
  double dbwdh = 0.0;
  double dbeffdh = 0.0;
  double dhfdh = 0.0;
  double dflcovdh = 0.0;
  double dwcovdh = 0.0;

  if (parameterID == 1) // d
    dddh  = 1.0;
  else if (parameterID == 2) // bw 
    dbwdh = 1.0;
  else if (parameterID == 3) // beff
    dbeffdh = 1.0;
  else if (parameterID == 4) // hf
    dhfdh =  1.0;
  else if (parameterID == 7) // flcover
    dflcovdh =  1.0;
  else if (parameterID == 8) // wcover
    dwcovdh =  1.0;
  else {
    for (int i = 0; i < nFibers; i++)
      dyidh[i] = 0.0;
  }

  if (parameterID == 1 || parameterID == 2 || parameterID == 3 || parameterID == 4 || parameterID == 7 || parameterID == 8  ) {

    double A = (d-hf)*(d-hf)*bw/2 + beff*hf*(d-hf/2);
    double B = beff*hf + (d-hf)*bw;
    double dAdh = 1/2*dbwdh*(d-hf)*(d-hf)+ 2*bw*(d-hf)*(dddh-dhfdh)+dbeffdh*hf*(d-hf/2)+ beff*dhfdh*(d-hf/2)+beff*hf*(dddh-dhfdh/2);
    double dBdh = dbeffdh*hf+dhfdh*beff+dbwdh*(d-hf)+bw*(dddh-dhfdh);
    double dYgdh = dAdh/B-A/B/B*dBdh;

    int loc = 0;
    
    // derive Location of core web fibers
    double dyIncrwdh  = (dddh-dhfdh-dwcovdh)/Nwcore;
    double dyStartwdh = -dYgdh + dwcovdh + 0.5 * dyIncrwdh;
    for (int i = 0; i < Nwcore; i++)
      dyidh[loc++] = dyStartwdh + dyIncrwdh*i;
    
    // derive Location of core flange fibers
    double dyIncrfldh = (dhfdh-dflcovdh) / Nflcore;
    double dyStartfldh = dddh - dYgdh -  dhfdh + 0.5* dyIncrfldh;
    for (int i = 0; i < Nflcore; i++)
      dyidh[loc++] = dyStartfldh + dyIncrfldh*i;
    
    // derive Location of cover web fibers
    double dyIncrwcovdh = dwcovdh / Nwcover;
    double dyStartwcovdh = -dYgdh + 0.5 *dyIncrwcovdh;
    for (int i = 0; i < Nwcover; i++)
      dyidh[loc++] = dyStartwcovdh + dyIncrwcovdh*i;
    
    // derive Location of cover flange fibers
    double dyIncrflcovdh = dflcovdh /Nflcover;
    double dyStartflcovdh = dddh - dYgdh - dflcovdh + 0.5*dyIncrflcovdh;
    for (int i = 0; i < Nflcover; i++)
      dyidh[loc++] = dyStartflcovdh + dyIncrflcovdh*i;
    
    // derive Location of Steel fibers
    double dyTopdh = dddh - dYgdh - dflcovdh;
    for (int i = 0; i < NsteelTop; i++)
      dyidh[loc++] = dyTopdh;
    
    double dyBotdh = -dYgdh + dwcovdh;
    for (int i = 0; i < NsteelBottom; i++)
      dyidh[loc++] = dyBotdh;
  }
  
  if (dzidh != 0) {
    for (int i = 0; i < nFibers; i++)
      dzidh[i] = 0.0;
  }

  return;
}

void
RCTBeamSectionIntegration::getWeightsDeriv(int nFibers, double *dwtsdh)
{
  double dddh  = 0.0;
  double dhfdh = 0.0;
  double dbwdh = 0.0;
  double dbeffdh = 0.0;
  double dflcovdh = 0.0;
  double dwcovdh = 0.0;
  double dAtopdh = 0.0;
  double dAbottomdh = 0.0;
  
  if (parameterID == 1)  // d --> 
    dddh  = 1.0;
  else if (parameterID == 2) // bw -->
    dbwdh = 1.0;
  else if (parameterID == 3) // beff
    dbeffdh = 1.0;
  else if (parameterID == 4) // hf
    dhfdh = 1.0;
  else if (parameterID == 5) // Atop
    dAtopdh = 1.0;
  else if (parameterID == 6) // Abottom
    dAbottomdh = 1.0;
  else if (parameterID == 7) // flcover
    dflcovdh = 1.0;
  else if (parameterID == 8) // wcover
    dwcovdh =  1.0;
  else {
    for (int i = 0; i < nFibers; i++)
      dwtsdh[i] = 0.0;
  }

  if (parameterID >= 1 && parameterID <= 8) {

    int loc = 0;
    
    // derive Weight of core web fibers
    double dAwcoredh = dbwdh*(d-hf-wcov)/Nwcore+bw*(dddh-dhfdh-dwcovdh)/Nwcore;
    for (int i = 0; i < Nwcore; i++) {
      dwtsdh[loc++] = dAwcoredh;
    }
    
    // derive Weight of core flange fibers
    double dAflcoredh = (dhfdh-dflcovdh)*beff/Nflcore+(hf-flcov)*dbeffdh/Nflcore;
    for (int i = 0; i < Nflcore; i++) {
      dwtsdh[loc++] = dAflcoredh;
    }
    
    // derive Weight of cover web fibers
    double dAwcoverdh = dbwdh * wcov/Nwcover + bw*dwcovdh/Nwcover;
    for (int i = 0; i < Nwcover; i++) {
      dwtsdh[loc++] = dAwcoverdh;
    }
    
    // derive Weight of cover flange fibers
    double dAflcoverdh = dbeffdh*flcov/Nflcover + dflcovdh*beff/Nflcover;
    for (int i = 0; i < Nflcover; i++) {
      dwtsdh[loc++] = dAflcoverdh;
    }
    
    // derive Weight of Steel fibers
    for (int i = 0; i < NsteelTop; i++)
      dwtsdh[loc++] = dAtopdh;

    for (int i = 0; i < NsteelBottom; i++)
      dwtsdh[loc++] = dAbottomdh;    
  }

  return;
}

void
RCTBeamSectionIntegration::Print(OPS_Stream &s, int flag)
{
  s << "RCT" << endln;
  s << " d = "  << d;
  s << " bw = " << bw; 
  s << " beff = " << beff;
  s << " hf = " << hf; 
  s << " Atop = " << Atop;
  s << " Abottom = " << Abottom;
  s << " flcov = " << flcov;
  s << " wcov = " << wcov << endln;
  s << " Nflcover = " << Nflcover;
  s << " Nwcover = " << Nwcover;
  s << " Nflcore = " << Nflcore;
  s << " Nwcore = " << Nwcore;
  s << " NsteelTop = " << NsteelTop;
  s << " NsteelBottom = " << NsteelBottom << endln;

  return;
}

int
RCTBeamSectionIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(14);

  data(0) = d;
  data(1) = bw;
  data(2) = beff;
  data(3) = hf;
  data(4) = Atop;
  data(5) = Abottom;
  data(6) = flcov;
  data(7) = wcov;
  data(8) = Nflcover;
  data(9) = Nwcover;
  data(10) = Nflcore;
  data(11) = Nwcore;
  data(12) = NsteelTop;
  data(13) = NsteelBottom;


  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "RCTBeamSectionIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
RCTBeamSectionIntegration::recvSelf(int cTag, Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
{
  static Vector data(14);
  
  int dbTag = this->getDbTag();
  
  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "RCTBeamSectionIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  d	       = data(0);
  bw	       = data(1);
  beff         = data(2);
  hf	       = data(3);
  Atop         = data(4);
  Abottom      = data(5);
  flcov        = data(6);
  wcov         = data(7);
  Nflcover     = (int)data(8);
  Nwcover      = (int)data(9);
  Nflcore      = (int)data(10);
  Nwcore       = (int)data(11);
  NsteelTop    = (int)data(12);
  NsteelBottom = (int)data(13);

  return 0;
}
