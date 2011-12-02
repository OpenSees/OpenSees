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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/integration/WideFlangeSectionIntegration.h,v $

#ifndef WideFlangeSectionIntegration_h
#define WideFlangeSectionIntegration_h

#define SECTION_INTEGRATION_TAG_WideFlange 1976

#include <SectionIntegration.h>

class UniaxialMaterial;

class WideFlangeSectionIntegration : public SectionIntegration
{
 public:
  WideFlangeSectionIntegration(double d, double tw, double bf, double tf,
			       int Nfdw, int Nftf);
  WideFlangeSectionIntegration();
  ~WideFlangeSectionIntegration();

  int getNumFibers(void);

  void getFiberLocations(int nFibers, double *xi);
  void getFiberWeights(int nFibers, double *wt);

  SectionIntegration *getCopy(void);

  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);

  void getLocationsDeriv(int nFibers, double *dptsdh);
  void getWeightsDeriv(int nFibers, double *dwtsdh);

  void Print(OPS_Stream &s, int flag = 0);

  int arrangeFibers(UniaxialMaterial **theMaterials,
		    UniaxialMaterial *theSteel);

 private:
  double d;
  double tw;
  double bf;
  double tf;

  int Nfdw;
  int Nftf;

  int parameterID;
};

#endif
