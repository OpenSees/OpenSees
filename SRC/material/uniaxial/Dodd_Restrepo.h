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

 // Written: fmk
 // Created: 10/2011
 //
 // Description: This file contains the class definition for 
 // Dodd_Restrepo. Dodd_Restrepo is a c++ wrapper to the steel material
 // found in DoddRestrepo.f that was written by L.L. Dodd and J. Restrepo  1991

 #ifndef Dodd_Restrepo_h
 #define Dodd_Restrepo_h

 #include <UniaxialMaterial.h>

 class Dodd_Restrepo : public UniaxialMaterial
 {
 public:
   Dodd_Restrepo(int tag, 
		 double Fy, 
		 double Fsu, 
		 double ESH, 
		 double ESU,
		 double Youngs, 
		 double ESHI, 
		 double FSHI, 
		 double OmegaFac,
		 double Conv);
   
   int setTrialStrain(double strain, double strainRate = 0.0);
   int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0);
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
   
   void Print(OPS_Stream &s, int flag = 0);
   
 protected:

 private:

   double FSH;
   double ESU;
   //   int invokeSubroutine(int ist);
   
   // input parameters
   double Fy       ; // Yield Stress
   double Fsu      ; // Engineering coordinate "ultimate" stress
   double ESH      ; // Engineering coordinate strain hardening strain
   double Youngs   ; // Youngs modulus
   double ESHI     ; // Intermediate strain hardening curve engineering strain
   double FSHI     ; // Intermediate strain hardening curve engineering stress
   double Conv     ; // strain factor (1 for strain 1.E6 for microstrain)
   double OmegaFac ; // Multiplication factor for Omega
   
   // addiional function input
   double Es       ; // Engineering Strain
   double EpsLast  ; // Natural strain at last increment
   double FpsLast  ; // True stress at last increment
   double YpTanLast; // Tangential modulus at last increment
   double EpsOld   ; // Natural strain at second to last increment
   double Epy      ; // The yield strain, Fy/Youngs (a positive value)
   double EpSH     ; // Natural coordinate strain hardening strain
   double Epsu     ; // Natural coordinate "ultimate" strain
   double Fpsu     ; // True coordinate "ultimate" stress (slope at ultimate)
   double SHPower  ; // Exponent which governs the strain-hardening curve
   
   // changeable
   double Epr[2]   ; // Reversal strain (1=tension, 2=compression)
   double Fpr[2]   ; // Reversal stress (1=tension, 2=compression)
   double Epa[2]   ; // Strain at end of linear branch (1=tension, 2=compression)
   double Fpa[2]   ; // Stress at end of linear branch (1=tension, 2=compression)
   double Epo[2]   ; // Maximum "natural" shift (1=compression, 2=tension)
   double EpoMax   ; // The maximum magnitude of Epo(1) and Epo(2)
   double EpsuSh[2]; // Shifted "ultimate" strain value (1=tension, 2=compression)
   double YoungsUn ; // Unloading modulus   // Invokes the FORTRAN subroutine
   double Power[2] ; // Exponent in normalised Bauschinger eq. (1=tens., 2=comp.)
   int BFlag[2];
   int  LMR;
   double EprM[2]  ; // Major reversal strain (1=tension, 2=compression)
   double FprM[2]  ; // Major reversal stress (1=tension, 2=compression)
   double EpaM[2]  ; // Major reversal Epa
   double FpaM[2]  ; // Major reversal Fpa
   double YpTanM[2]; // Tangential modulus at major reversals (1=tens, 2=comp)
   double PowerM[2]; // Major reversal Power
   
   
   // output
   double Eps      ; // Natural strain
   double Fps      ; // True coordinates stress
   double Fs       ; // Engineering Stress
   double YpTan    ; // True coordinates tangential modulus
   double YTan     ; // Tangential modulus
   
   double tStrain;
   double tTangent;
   double tStress;
   double cStrain;
   double cTangent;
   double cStress;
   
   double cEpr[2]   ; // Reversal strain (1=tension, 2=compression)
   double cFpr[2]   ; // Reversal stress (1=tension, 2=compression)
   double cEpa[2]   ; // Strain at end of linear branch (1=tension, 2=compression)
   double cFpa[2]   ; // Stress at end of linear branch (1=tension, 2=compression)
   double cEpo[2]   ; // Maximum "natural" shift (1=compression, 2=tension)
   double cEpoMax   ; // the maximum magnitude of Epo(1) and Epo(2)
   double cEpsuSh[2]; // Shifted "ultimate" strain value (1=tension, 2=compression)
   double cYoungsUn ; // Unloading modulus   // Invokes che FORTRAN subroutine
   double cPower[2] ; // Exponent in normalised Bauschinger eq. (1=tens., 2=comp.)
   int    cBFlag[2];
   int    cLMR;
   double cEprM[2]  ; // Major reversal strain (1=tension, 2=compression)
   double cFprM[2]  ; // Major reversal stress (1=tension, 2=compression)
   double cEpaM[2]  ; // Major reversal Epa
   double cFpaM[2]  ; // Major reversal Fpa
   double cYpTanM[2]; // tangential modulus at major reversals (1=tens, 2=comp)
   double cPowerM[2]; // Major reversal Power

   int myTag;
};

#endif
