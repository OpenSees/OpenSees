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
                                                                        
// $Revison
// $Date
// $Source

// written by: Davide Fugazza

// Description: This file contains the class definition for SMAMaterial. 

#ifndef SMAMaterial_h
#define SMAMaterial_h

#include <UniaxialMaterial.h>

class SMAMaterial : public UniaxialMaterial
{
 public:
  SMAMaterial(int tag, double E, double eps_L, double sig_AS_s, double sig_AS_f, double sig_SA_s, double sig_SA_f);
  SMAMaterial();    
  ~SMAMaterial();
  
  int setTrialStrain(double strain, double strainRate = 0.0); 
  double getStrain(void);
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void);
  
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

  // Material parameters
  double E;        // elastic modulus for austenite and martensite 
  double eps_L;    // plateau length (ex. 0.08 means 8%)                              
  double sig_AS_s; // forward transformation stress start          
  double sig_AS_f; // forward transformation stress finish
  double sig_SA_s; // reverse transformation stress start
  double sig_SA_f; // reverse transformation stress finish
  
  // Committed history variables
  double Cstrain;  // strain              at time t_n (previous time step)
  double Cstress;  // stress              at time t_n (previous time step)
  double Ccsi;     // martensite fraction at time t_n (previous time step)
  
  // Trial history variables
  double Tcsi;     // trial martensite fraction at time t (current time step)
  
  // Trial state variables
  double Tstrain;  // trial strain at time t (current time step)
  double Tstress;  // trial stress at time t (current time step)
  double Ttangent; // tangent      at time t (current time step)
};

#endif

