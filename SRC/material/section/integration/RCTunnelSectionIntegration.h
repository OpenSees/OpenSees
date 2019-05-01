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

// $Revision$
// $Date$
// $Source$

#ifndef RCTunnelSectionIntegration_h
#define RCTunnelSectionIntegration_h

#include <SectionIntegration.h>

class UniaxialMaterial;

class RCTunnelSectionIntegration : public SectionIntegration
{
 public:
  RCTunnelSectionIntegration(double d, double h, double Asin, double Asout, 
			     double coverin, double coverout,
			     int Nrings, int Nwedges, int Nbarsin, int Nbarsout);
  RCTunnelSectionIntegration();
  ~RCTunnelSectionIntegration();

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

  int arrangeFibers(UniaxialMaterial **theMaterials,
		    UniaxialMaterial *theConcrete,
		    UniaxialMaterial *theSteel);

 private:
  double d; // outer diameter
  double h; // lining thickness
  double Asinner; // bar area inner ring
  double Asouter; // bar area outer ring
  double coverinner; // cover to inner ring
  double coverouter; // cover to outer ring

  int Nrings;
  int Nwedges;
  int Nbarsinner;
  int Nbarsouter;

  int parameterID;
};

#endif
