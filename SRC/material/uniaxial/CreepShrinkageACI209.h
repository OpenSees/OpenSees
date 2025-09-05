/* ****************************************************************** **
** OpenSees - Open System for Earthquake Engineering Simulation       **
** Pacific Earthquake Engineering Research Center                     **
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
** Frank McKenna (fmckenna@ce.berkeley.edu)                         **
** Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
** Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
   
 //----------------------------------------------------------------------------------------------------------------------------
 // Developed by:                                                                
 // Javad Esmaeelpour (jesmaeel@tennessee.edu)    
 // Mark D. Denavit   (mdenavit@utk.edu)           
 // Michael H. Scott  (michael.scott@oregonstate.edu)
 //
 // Based on TDConcrete implementations by:
 // Adam M. Knaack (adam.knaack@schaefer-inc.com) 
 // Schaefer-Inc, Cincinnati, Ohio, USA
 // Nikola D. Tosic (ntosic@imk.grf.bg.ac.rs)
 // Department for Materials and Structure, Faculty of Civil Engineering, University of Belgrade, Serbia
 // Yahya C. Kurama (ykurama@nd.edu)
 // Department of Civil and Environmental Engineering and Earth Sciences, College of Engineering, University of Notre Dame, Notre Dame, Indiana, USA
 //----------------------------------------------------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------------------------------------------------
 // Description: This file contains the source code of CreepShrinkageACI209. 
 // CreepShrinkageACI209 is a wrapper that imposes creep and shrinkage evoluation equations
 // to any uniaxialMaterial.
 //----------------------------------------------------------------------------------------------------------------------------

#ifndef CreepShrinkageACI209_h
#define CreepShrinkageACI209_h

#include <UniaxialMaterial.h>
#include <Domain.h>

class CreepShrinkageACI209 : public UniaxialMaterial
{
public:
  // -- Core UniaxialMaterial Interface
  CreepShrinkageACI209(int tag, double _fc, double _fcu, double _epscu, double _Ec, double _age, double _epsshu, double _epssha, double _tcr, double _epscru, double _epscra, double _epscrd, double _tcast);
  CreepShrinkageACI209(int tag, UniaxialMaterial &matl, double _age, double _epsshu, double _epssha, double _tcr, double _epscru, double _epscra, double _epscrd, double _tcast);
  CreepShrinkageACI209(void);
  virtual ~CreepShrinkageACI209();
  
  const char *getClassType(void) const {return "CreepShrinkageACI209";};    
  double getInitialTangent(void);
  UniaxialMaterial *getCopy(void);
  
  int setTrialStrain(double strain, double strainRate = 0.0); 
  double getStrain(void);
  double getStress(void);
  double getTangent(void);

  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);        
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
  int getVariable(const char *variable, Information &);
  
  // -- Functions for Recorders
  Response *setResponse(const char **argv, int argc,OPS_Stream &theOutput);
  int getResponse(int responseID, Information &matInfo);
  
protected:
  
private:
  // -- Private Functions (Implementation Details)
  double setCreepStrain(double time);    // creep strain at global time
  double getCurrentTime(void);           // Current analysis time from the active OpenSees Domain.
  double setShrink(double time);         // Shrinkage strain at global time
  void expandArrays(void);               // Ensures capacity for the stress/time history arrays.
   

  // -- Member Variables
  UniaxialMaterial *wrappedMaterial;
  
  // matpar : Concrete FIXED PROPERTIES

  double tcr;                        // Tcr, creep age at loading used for phiu
  double Ec;                         // Concrete modulus at loading
  double age;                        // tD, analysis time when drying starts
  double epsshu;                     // (eps_sh)u, ultimate shrinkage strain (ACI 209R Eq. 2-7)
  double epssha;                     // f, shrinkage time-shape parameter (ACI 209R Eq. 2-7)
  double epscra;                     // psi, creep time-shape parameter (ACI 209R Eq. 2-6)
  double epscru;                     // phiu, ultimate creep coefficient (ACI 209R Eq. 2-6)
  double epscrd;                     // d, creep time-shape parameter (ACI 209R Eq. 2-6)
  double tcast;                      // Casting time in global analysis clock
  
  
  double trialStress;                // Current trial stress from wrapped material
  double trialTangent;               // Current trial tangent from wrapped material
  double trialTotalStrain;           // Total applied strain (input)
  double trialCreepStrain;           // Computed creep strain component
  double trialShrinkageStrain;       // Computed shrinkage strain component 
  double trialMechanicalStrain;      // Mechanical strain sent to wrapped material
  
  double committedStress;		         // Last committed stress
  double committedTangent;           // Last committed tangent
  double committedCreepStress;       // Last committed creep stress
  double committedCreepStrain;       // Last committed creep strain
  double committedShrinkageStrain;   // Last committed shrinkage strain
  double committedMechanicalStrain;  // Last committed mechanical strain
  double committedTotalStrain;       // Last committed total strain

  int    historyPointCount;          // Number of points in stress-time history
  int    iterationInStep;            // Iteration counter to avoid recomputing within same step 

  enum{startSize=500, growSize=200};
  int maxSize;                       // Current allocated size of history arrays
  
  double *DSIG_i;                    // Stress increment history array [historyPointCount]
  double *TIME_i;                    // Time history array [historyPointCount]
};


#endif