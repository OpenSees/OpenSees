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
                                                                        
// Concrete Stress-Strain Model Including Unloading/Reloading Model
// Proposed by Sakai and Kawashima (Envelope Curve: Mander Model)
// Coded by J. Sakai    2002. 2.14

// modified for OpenSees: fmk 6/14

#ifndef ConcreteSakaiKawashima_h
#define ConcreteSakaiKawashima_h

#include <UniaxialMaterial.h>

class ConcreteSakaiKawashima : public UniaxialMaterial
{
  public:
  ConcreteSakaiKawashima(int tag, double YMx, double Sigcc, double EPScc);
  ConcreteSakaiKawashima(void);
  ~ConcreteSakaiKawashima();

  const char *getClassType(void) const {return "ConcreteSakaiKawashima";};

  int setTrialStrain(double strain, double strainRate = 0.0); 
  double getStrain(void);      
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void) {return YMc;}

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
  // input
  double YMc;    // initial Stiffness
  double Sigcc;  // peak stress
  double EPScc;  // strain at peak stress

  // state variables
  double tTangent;
  double tStrain;
  double tStress;
  double cTangent;
  double cStrain;
  double cStress;


  double DE0;    // Delta Epsilon of the Previous Step     
  double Sigule; // Stress at the Unloading Point on Envelope
  double EPSule; // Strain at the Unloading Point on Envelope 
  double Sigul;  // Stress at the Unloading Point
  double EPSul;  // Strain at the Unloading Point                       
  double EPSpl;  // Plastic Strain
  double Suln;   // Stress at the Unloading Point after Un/Reloading
  double YMrl;   // Slope of Reloding Branch
  double YMtan;  // Tangential Slope of the Previous Step
  double Sigrl;  // Stress at the Reloading Point 
  double EPSrl;  // Strain at the Reloading Point
  double EPSpl0; // EPSpl of Previous Unloading 
  double Suln0;  // Suln of Previous Reloading 
  double GamRL;  // Partial Reloading Ratio GammaRL
  int Jcon;      // Index of Current Status 
  int Ncyc;
  int Jcon0;
  int Ncyc0;

  // committed state variables
  double cDE0;    // Delta Epsilon of the Previous Step     
  double cSigule; // Stress at the Unloading Point on Envelope
  double cEPSule; // Strain at the Unloading Point on Envelope 
  double cSigul;  // Stress at the Unloading Point
  double cEPSul;  // Strain at the Unloading Point                       
  double cEPSpl;  // Plastic Strain
  double cSuln;   // Stress at the Unloading Point after Un/Reloading
  double cYMrl;   // Slope of Reloding Branch
  double cYMtan;  // Tangential Slope of the Previous Step
  double cSigrl;  // Stress at the Reloading Point 
  double cEPSrl;  // Strain at the Reloading Point
  double cEPSpl0; // EPSpl of Previous Unloading 
  double cSuln0;  // Suln of Previous Reloading 
  double cGamRL;  // Partial Reloading Ratio GammaRL
  int cJcon;      // Index of Current Status 
  int cNcyc;
  int cJcon0;
  int cNcyc0;
};


#endif

