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
// $Date: 2008-04-14 21:20:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete04.h,v $
                                                                        
                                                                        
#ifndef Concrete04_h
#define Concrete04_h

// File: Concrete04.h
//
// Written: N.Mitra (nmitra@u.washington.edu)
// Created: 09/04
// Revision: A
//
// Description: This file contains the class definition for 
// Concrete04.h 
//   - No tension 
//   - Linear unloading/reloading
//
// What: "@(#) Concrete04.h, revA"
// Revision: 1. Adding in Exponential tension part (05-16-05)

#include <UniaxialMaterial.h>

class Concrete04 : public UniaxialMaterial
{
 public:
//  Concrete04 (int tag, double fpc, double eco, double ecu, double Ec0, double fct);
  Concrete04 (int tag, double fpc, double eco, double ecu, double Ec0, double fct, double etu);  
  Concrete04 (int tag, double fpc, double eco, double ecu, double Ec0, double fct, double etu, double beta);
  Concrete04 (int tag, double fpc, double eco, double ecu, double Ec0);
  Concrete04 ();
  ~Concrete04();

  const char *getClassType(void) const {return "Concrete04";};
  
  int setTrialStrain(double strain, double strainRate = 0.0); 
  double getStrain(void);      
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void) {return Ec0;}
  
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);        
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
  // LOWES: add function for use with variable hinge lnegth model
  int getMaterialType(void);
  // LOWES: end
  
 protected:
  
 private:
  /*** Material Properties ***/
  double fpc;    // Compressive strength
  double epsc0;  // Strain at compressive strength
  double epscu;  // Strain at crushing strength
  double Ec0;    // initial tangent
  double fct;   // Concrete tensile strength
  double etu;   // ultimate tensile strain              
  double beta;  // exponential curve parameter, residual stress (as a factor of ft)
  // at etu. 
  
  /*** CONVERGED History Variables ***/
  double CminStrain;   // Smallest previous concrete strain (compression)
  double CmaxStrain;  
  double CunloadSlope; // Unloading (reloading) slope from CminStrain
  double CendStrain;   // Strain at the end of unloading from CminStrain
  double CcompStrain;   // strain value at which the compression unloading intersects the   
  // zero stress value or the strain value at which tensile reloading starts.                                  
  double CUtenStress;      // tensile stress value at which unloading begins
  double CUtenSlope;      // unloading tensile slope value
  
  /*** CONVERGED State Variables ***/
  double Cstrain;
  double Cstress;   
  double Ctangent;	// Don't need Ctangent other than for revert and sendSelf/recvSelf
  // Storing it is better than recomputing it!!!
  
  double TminStrain;
  /*** TRIAL History Variables ***/      
  double TmaxStrain;
  double TunloadSlope;
  double TendStrain;
  double TcompStrain;
  double TUtenStress;
  double TUtenSlope;
  
  /*** TRIAL State Variables ***/
  double Tstrain;
  double Tstress;
  double Ttangent; // Not really a state variable, but declared here
  // for convenience
  
  void CompReload(void);
  void CompEnvelope(void);
  void setCompUnloadEnv(void);
  void TensReload(void);
  void TensEnvelope(void);
  void setTenUnload(void);
};


#endif
