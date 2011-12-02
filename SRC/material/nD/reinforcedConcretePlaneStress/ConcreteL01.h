#ifndef ConcreteL01_h
#define ConcreteL01_h

// File ConcreteL01.h
// Hsu and Mansour's Model
// Written: alaskar
// Created: 2007.8
// For Detailed explanation of the model, please refer to the book
// entitled "Unified Theory of Concrete Structures,"
// by Thomas T.C. Hsu and Y.L. Mo, John Wiley & Sons, April 2010.

#include <UniaxialMaterial.h>

class ConcreteL01 : public UniaxialMaterial
{
 public:
  ConcreteL01(int tag, double fpc, double eco);
  ConcreteL01();
  ~ConcreteL01();
  
  int setTrialStrain(double strain, double strainRate = 0.0); 	
  int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
  
  //  int setTrialStrain(double x, double k, double Dfactor, double BETA, double EPSLONTP, double strain, double strainRate = 0.0);
  // int setTrial (double x, double k, double Dfactor, double BETA, double EPSLONTP, double strain, double &stress, double &tangent, double strainRate=0.0);

  //beta and epslonTP are required arguments when to calculate softening effect zeta
  // D: damage factor for strength
  // x, k are required for the improvement by the normal stresses
  
  double getStrain(void);              
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void) {return 1.7*fpc/epsc0;};
  double getSecant(void);
  double getZeta(void);
  double getPD(void); // Get partial differentiation of stress to epslonTP (strain in perpendicular direction)
  
  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);   
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    

  Response *setResponse (const char **argv, int argc, 
			 OPS_Stream &theOutputStream);
  int getResponse (int responseID, Information &matInformation);      

  void Print(OPS_Stream &s, int flag =0);
  
 protected:
  
 private:
  
  /*** Material Properties ***/        
  double fpc;      // Compressive strength
  double epsc0;    // Strain at compressive strength
  double zeta;     // Softening effect
  double beta;     // Parameter needed for calculating zeta
  double fbeta;    // function of beta
  double Wp;		 // prestressing factor
  double epslonTP; // Strain in the perpendicular direction, needed to get the zeta
  //double epslonCI; // initial prestressing strain
  //double sigmaCI;  // initial prestress
  double D;        // Damage factor for strength, get from parameter
  double X;        // for normal stresses 
  double K;        // for normal stresses
  
  /*** History Variables ***/
  int CloadingState; // Flag for loading state
  // 1 = ascending branch of envelope in compression
  // 2 = descending branch of envelope in compression
  // 3 = ascending branch of envelope in tension
  // 4 = descending branch of envolope in tension
  // 5 = reloading at compression envelope 1 or 2
  // 6 = unloading from tenison to compression zone ( from 4 to 1 or 2 )
  // 7 = reloading from compression to tension (reload from 5 approach to tension zone)

  
  int reloadPath;   // Flag for reversing state
  // = 0 initially
  // = 1 if reverse from ascending compression branch
  // = 2 if reverse from descending compression branch
  
  double reverseFromOneStrain; // Strain and stress of the reversed point from path 1
  double reverseFromOneStress;
  double reverseFromTwoStrain; // Strain and stress of the reversed point from path 2
  double reverseFromTwoStress;
  double reverseFromFourStrain; // Strain and stress of the reversed point from path 4
  double reverseFromFourStress;
  
  /*** Current Key Points According To State Variables ***/
  double interFiveSevenStrain;  // Strain of start point of path 7 (end point of path 5)
  
  double approachFiveToComStrain;    // Intersection point of path 5 to compressive envelope
  double approachSixToComStrain;     // Intersection point of path 6 to compressive envelope
  
  
  /*** CONVERGED State Variables ***/    
  double Cstrain;
  double Cstress;
  double Ctangent;    
  
  /*** TRIAL History Variables ***/
  int TloadingState;	
  
  /*** TRIAL State Variables ***/
  double Tstrain;
  double Tstress;
  double Ttangent; // Not really a state variable, but declared here
  // for convenience
  
  
  // Calculates the trial state variables based on the trial strain
  void determineTrialState ( double dStrain);
  
  //Envelope curve in beginning
  void envelope( );
  
  void pathFive( );
  void pathSix ( );
  void pathSeven ();
  void getApproachFiveToComStrain ( );
  void getApproachSixToComStrain ( );
};

#endif
