                                                                
                                                                       
#ifndef TendonL01_h
#define TendonL01_h

//TendonL01.h, Hsu and Wang's Model
//Written:  ALaskar
//Created:  2007.8
//Description: This file contains the class definition for TendonL01.h
// For Detailed explanation of the model, please refer to the book
// entitled "Unified Theory of Concrete Structures,"
// by Thomas T.C. Hsu and Y.L. Mo, John Wiley & Sons, April 2010.
#include <UniaxialMaterial.h>

// Default values for reloading/unloading path parameters AC, RC
#define TENDON_01_DEFAULT_AC       1.9
#define TENDON_01_DEFAULT_RC       10.0
#define LOOP_NUM_LIMIT               30
const int SIZE1 = LOOP_NUM_LIMIT; //limit of array number

class TendonL01 : public UniaxialMaterial
{
  public:
    TendonL01(int tag, double fy, double E0, double fpu, double rou, double epsp,
	      double ac = TENDON_01_DEFAULT_AC,
	      double rc = TENDON_01_DEFAULT_RC);
    TendonL01();
    ~TendonL01();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return Eps;};
    double getSecant(void);

    double getCommittedStrain(void);
    
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
    double fpy;		// yield strength of bare steel
    double Eps;		// initial stiffness of steel
    double fpu;		// compressive strength of concrete	
    double rou;		// steel ration
    double epsp;		// prestressinf strain
    double ac;		// parameter of reverse loop
    double rc;		// parameter of reverse loop
    
    /*** CONVERGED History Variables ***/
    double CminStrain;  // Minimum strain in compression
    double CmaxStrain;  // Maximum strain in tension
    int CloadingState; // Flag for loading state
    // 1 = initial envelope
    // 2 = tension envelope after yield
    // 3 = compression envelope after yield
    // 4 = reverse from tension envelope
    // 5 = reverse from compression envelope
    int CloopPathState; // Flag for loop state
    // 1, 2 and 3 for the loop going down
    // 1 : the first path when stress >0
    // 2 : the second path when stress <0 to -0.65fy
    // 3 : the third path when stress <-0.65fy to the end point of the loop
    // 4, 5 and 6 for the loop going up
    // 4 : the first path when stress <0
    // 5 : the second path when stress >0 to 0.65fy
    // 6 : the third path when stress >0.65fy to the end point of the loop
    
    double reverseFromTenEnvelopeStrain;  // Strain of reversed point from tensile envelope
    double reverseFromTenEnvelopeStress;  // Stress of reversed point from tensile envelope
    double approachToComEnvelopeStrain;   // Strain of point when unloading path approaches to compression envelope
    double approachToComEnvelopeStress;   // Stress of point when unloading path approaches to compression envelope
    
    double reverseFromComEnvelopeStrain;  // Strain of reversed point from compressive envelope
    double reverseFromComEnvelopeStress;  // Stress of reversed point from compressive envelope
    double approachToTenEnvelopeStrain;   // Strain of point when reloading path approaches to tensile envelope
    double approachToTenEnvelopeStress;   // Stress of point when reloading path approaches to tensile envelope
    
    // History converged variables for hysteristic loop
    double CreverseTopStrain[SIZE1];    // array storing the strain and stress of top point of reverse loop
    double CreverseTopStress[SIZE1];
    double CreverseBottomStrain[SIZE1]; // array storing the strain and stress of last point of reverse loop
    double CreverseBottomStress[SIZE1];
    int    CreverseTopNum;             // Num of reverse points reversed from up
    int    CreverseBottomNum;          // Num of reverse points reversed from down
    
    // Trial variables for hysteristic loop
    double TreverseTopStrain[SIZE1];    // array storing the strain and stress of top point of reverse loop
    double TreverseTopStress[SIZE1];
    double TreverseBottomStrain[SIZE1]; // array storing the strain and stress of last point of reverse loop
    double TreverseBottomStress[SIZE1];
    int    TreverseTopNum;             // Num of reverse points reversed from up
    int    TreverseBottomNum;          // Num of reverse points reversed from down
    
    double downPathPointOneStrain;    // key points at down path of loop
    double downPathPointOneStress;    // set to be zero at default value
    double downPathPointTwoStrain;
    double downPathPointTwoStress;
    
    double upPathPointOneStrain;      // key points at up path of loop
    double upPathPointOneStress;      // set to be zero at default value
    double upPathPointTwoStrain;
    double upPathPointTwoStress;
    
    
    /*** CONVERGED State Variables ***/    
    double Cstrain;
    double Cstress;
    double Ctangent;    
    
    /*** TRIAL History Variables ***/
    double TminStrain;  // Minimum strain in compression
    double TmaxStrain;  // Maximum strain in tension
    int TloadingState;  // Flag for loading state, definition is same as CloadingState
    int TloopPathState; // Flag for loop state, definition is same as CloopPathState
    
    
    /*** TRIAL State Variables ***/
    double Tstrain;
    double Tstress;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience
    
    
    void determineTrialState (double dStrain);  
    // Calculates the trial state variables based on the trial strain
    
    void initialEnvelope( ); 
    void tensionEnvelope( );
    void compressionEnvelope( );
    
    void determineTrialLoop(double dStrain); // Calacuate the trail state variables in the loop
    void determineDownPathPoint( ); //determine key points of down and up path of hysteristic loop
    void determineUpPathPoint( );
    
    void downPath( ); 
    void upPath( );
    
    void reverseFromTenEnvelope( );
    void reverseFromComEnvelope( );
    
    void reverseLoopSetZero( ); // Loop variables set to zero when merge into envelope
    
    double tt1; // for check
    double tt2;
    double ttStrain;
};

#endif
