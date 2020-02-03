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
                                                                        
// $Revision: 1.6 $
// $Date: 2007-05-07 21:30:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FatigueMaterial.h,v $
                                                      
// Written: Patxi
// Created: Aug 2003
//
// Description: This file contains the class definition for 
// FatigueMaterial.  FatigueMaterial wraps a UniaxialMaterial
// and imposes fatigue limits. More information about this material can
// be found in the doctoral dissertation of Patxi Uriz:
//
//   Uriz, Patxi, "Towards Earthquake Resistant Design of 
//      Concentrically Braced Steel Frames," Ph.D. Dissertation, 
//      Structural Engineering, Mechanics, and Materials, Civil 
//      and Envrironmental Engineering, University of California, 
//      Berkeley, December 2005
//



#ifndef FatigueMaterial_h
#define FatigueMaterial_h

#include <UniaxialMaterial.h>

class FatigueMaterial : public UniaxialMaterial
{
 public:
  // Default calibrated values from Ballio and Castiglioni Calibrations for 
  // European steel, wide flange sections.

  FatigueMaterial(int tag, UniaxialMaterial &material, 
		  double Dmax    =  1.0,
		  double E0      =  0.191,
		  double m       = -0.458,
		  double minStrain = -1.0e16,
		  double maxStrain =  1.0e16 );
  
  FatigueMaterial();
  ~FatigueMaterial();

  const char *getClassType(void) const {return "FatigueMaterial";};
  
  int setTrialStrain(double strain, double strainRate = 0.0); 
  double getStrain(void);          
  double getStrainRate(void);
  double getStress(void);
  double getTangent(void);
  double getDampTangent(void);
  double getInitialTangent(void) {return theMaterial->getInitialTangent();}
  
  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);        
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);

  Response *setResponse (const char **argv, int argc, OPS_Stream &s);
  int getResponse (int responseID, Information &matInformation);    
  bool hasFailed(void);  

  //by SAJalali
  virtual double getEnergy(void) { return energy; }
protected:
  
 private:
	 double energy, CStress; //SAJalali

  UniaxialMaterial *theMaterial;
  
  double DI; //Damage index
  double  X; //Range in consideration
  double  Y; //Previous Adjacent Range
  double  A; //Peak or valley 1
  double  B; //Peak or valley 2
  double  C; //Peak or valley 2
  double  D; //Peak or valley 4
  int   PCC; /*Previous Cycle counter flag if >1 then previous 'n' 
	       cycles did not flag a complete cycle */
  int   R1F; //Flag for first  peak count
  int   R2F; //Flag for second peak count
  double cSlope; //Current Slope
  double PS; //Previous slope
  double EP; //Previous Strain
  int    SF; /*Start Flag = 0 if very first strain, 
	       (i.e. when initializing)    = 1 otherwise */
  double DL; //Damage if current strain was last peak.
  
  double Dmax;
  double E0;
  double m;
  
  double minStrain;
  double maxStrain;
  
  bool Cfailed;
  double trialStrain;

  // added 6/9/2006
  // For recording strain ranges (SRXX) and Number of Cycles (NCXX)
  double SR1;  // Committed strain range at peak
  double NC1;  // Committed number of cycles at SR1 (i.e. 1.0 or 0.5)
  double SR2;  // Committed strain range 2 at PSEUDO peak - there are potentially two ranges
  double NC2;  // Committed number of cycles at SR2 2 (at PSEUDO peak) - there are potentially two ranges
  double SR3;  // Committed strain range 3 at PSEUDO peak - there are potentially two ranges
  double NC3;  // Committed number of cycles at SR2 3 (at PSEUDO peak) - there are potentially two ranges
  
};


#endif

