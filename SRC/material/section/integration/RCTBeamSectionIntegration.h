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
// $Date: 2007/01/25 19:53:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/integration/RCSectionIntegration.h,v $

#ifndef RCTBeamSectionIntegration_h
#define RCTBeamSectionIntegration_h

#include <SectionIntegration.h>

class NDMaterial;
class UniaxialMaterial;

class RCTBeamSectionIntegration : public SectionIntegration
{
 public:
  RCTBeamSectionIntegration(double d, double bw, double beff, double hf, double Atop, double Abottom, double flcov,
			    double wcov, int Nflcover, int Nwcover, int Nflcore, int Nwcore, int NsteelTop, int NsteelBottom);
  RCTBeamSectionIntegration();
  ~RCTBeamSectionIntegration();
  
  int getNumFibers(FiberType type = all);
  
  void getFiberLocations(int nFibers, double *yi, double *zi = 0);
  void getFiberWeights(int nFibers, double *wt);
  
  SectionIntegration *getCopy(void);
  
  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);
  
  void getLocationsDeriv(int nFibers, double *dyidh, double *dzidh = 0);
  void getWeightsDeriv(int nFibers, double *dwtdh);
  
  void Print(OPS_Stream &s, int flag = 0);
  
  int arrangeFibers(UniaxialMaterial **theUni, NDMaterial **theND,
		    NDMaterial *theCore,
		    NDMaterial *theCover,
		    UniaxialMaterial *theSteel);
  int arrangeFibers(UniaxialMaterial **theUni, 
		    UniaxialMaterial *theCore,
		    UniaxialMaterial *theCover,
		    UniaxialMaterial *theSteel);

 private:
  double d;
  double bw;
  double beff;
  double hf;
  double Atop;
  double Abottom;
  double flcov;
  double wcov;
  
  int Nflcover;
  int Nwcover;
  int Nflcore;
  int Nwcore;
  int NsteelTop;
  int NsteelBottom;
  
  int parameterID;
};

#endif
