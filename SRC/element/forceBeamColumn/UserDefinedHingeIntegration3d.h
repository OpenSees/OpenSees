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

// $Revision: 1.2 $
// $Date: 2003-03-15 00:09:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/UserDefinedHingeIntegration3d.h,v $

#ifndef UserDefinedHingeIntegration3d_h
#define UserDefinedHingeIntegration3d_h

#include <BeamIntegration.h>

#include <Vector.h>

class Channel;
class FEM_ObjectBroker;

class UserDefinedHingeIntegration3d : public BeamIntegration
{
 public:
  UserDefinedHingeIntegration3d(int npL, const Vector &ptL, const Vector &wtL,
				int npR, const Vector &ptR, const Vector &wtR,
				double E, double A, double Iz,
				double Iy, double G, double J);
  UserDefinedHingeIntegration3d();
  ~UserDefinedHingeIntegration3d();
  
  void getSectionLocations(int numSections, double L, double *xi);
  void getSectionWeights(int numSections, double L, double *wt);

  void addElasticDeformations(ElementalLoad *theLoad, double loadFactor,
			      double L, double *v0);
  int addElasticFlexibility(double L, Matrix &fe);

  BeamIntegration *getCopy(void);

  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  int setParameter(const char **argv, int argc, Information &info);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);

 private:
  Vector ptsL;
  Vector wtsL;
  Vector ptsR;
  Vector wtsR;

  double E;
  double A;
  double Iz;
  double Iy;
  double G;
  double J;
};

#endif
