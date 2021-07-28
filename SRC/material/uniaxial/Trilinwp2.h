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
                                                                        
// $Revision: 1.8 $
// $Date: 2008-12-18 23:40:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Trilinwp.h,v $

#ifndef Trilinwp2_h
#define Trilinwp2_h

#include <UniaxialMaterial.h>

class Trilinwp2 : public UniaxialMaterial
{
 public:
  Trilinwp2(int tag, 
		     double mom1p, double rot1p, double mom2p, double rot2p,
		     double mom3p, double rot3p,
		     double pinchX, double pinchY,
		     double damfc1 = 0.0, double damfc2 = 0.0,
		     double beta = 0.0, 
			double pt = 0.0, double pb = 0.0, double pc = 0.0, double mb = 0.0, int itype = 0);
  Trilinwp2();
  ~Trilinwp2();

  const char *getClassType(void) const {return "Trilinwp2";};
  
  int setTrialStrain(double strain, double strainRate = 0.0  );
  double getStrain(void);
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void) {return E1p;};
  
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
 protected:
  
 private:
  // Pinching parameters
  double pinchX;		// Deformation pinching
  double pinchY;		// Force pinching
  
  // Damage parameters
  double damfc1;		// Deformation
  double damfc2;		// Energy
  
  // Unloading parameter
  double beta;
  // Envelope parameters

  double pt;
  double pb;
  double pc;
  double mb;

  // Trial history variables
  double TrotMax;
  double TrotMin;
  double TrotPu;
  double TrotNu;
  double TenergyD;
  int TloadIndicator;
  
  // Trial state variables
  double Ttangent;
  double Tstress;
  double Tstrain;
  
  // Converged history variables
  double CrotMax;
  double CrotMin;
  double CrotPu;
  double CrotNu;
  double CenergyD;
  int CloadIndicator;
  
  // Converged state variables
  double Cstress;
  double Cstrain;
  
  // Backbone parameters
  double mom1p, rot1p;
  double mom2p, rot2p;
  double mom3p, rot3p;
  double mom1n, rot1n;
  double mom2n, rot2n;
  double mom3n, rot3n;
  
  double E1p, E1n;
  double E2p, E2n;
  double E3p, E3n;
  double Eup, Eun;

  double energyA;

  double mom1pi,mom2pi,mom3pi;
  double mom1ni, mom2ni, mom3ni;
  double rot1pi, rot2pi, rot3pi;
  double duct;
  int itype;

  void setEnvelope(void);
  
  double posEnvlpStress(double strain);
  double negEnvlpStress(double strain);
  
  double posEnvlpTangent(double strain);
  double negEnvlpTangent(double strain);
  
  double posEnvlpRotlim(double strain);
  double negEnvlpRotlim(double strain);
  
  void positiveIncrement(double dStrain);
  void negativeIncrement(double dStrain);
};

#endif
