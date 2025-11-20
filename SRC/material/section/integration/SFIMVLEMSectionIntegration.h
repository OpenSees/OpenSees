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
// $Date: 2010-09-13 21:31:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/integration/SFIMVLEMSectionIntegration.h,v $

#ifndef SFIMVLEMSectionIntegration_h
#define SFIMVLEMSectionIntegration_h

#include <SectionIntegration.h>

class NDMaterial;

class SFIMVLEMSectionIntegration : public SectionIntegration
{
 public:
  SFIMVLEMSectionIntegration(int Nlayers, double *h, double *b);
  SFIMVLEMSectionIntegration();
  ~SFIMVLEMSectionIntegration();

  int getNumFibers(FiberType type = all);

  void getFiberLocations(int nFibers, double *yi, double *zi = 0);
  void getFiberWeights(int nFibers, double *wt);

  const char* getClassType(void) const {return "SFIMVLEMSectionIntegration";}
  SectionIntegration *getCopy(void);

  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);

  void getLocationsDeriv(int nFibers, double *dyidh, double *dzidh = 0);
  void getWeightsDeriv(int nFibers, double *dwtdh);

  void Print(OPS_Stream &s, int flag = 0);

  int arrangeFibers(NDMaterial **theMaterials,
		    NDMaterial **theReinforcedConcrete);

 private:
  int Nlayers;
  
  double *h;
  double *b;

  int parameterID;
};

#endif
