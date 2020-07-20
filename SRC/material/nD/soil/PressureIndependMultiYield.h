// $Revision: 1.20 $
// $Date: 2009-01-16 19:40:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PressureIndependMultiYield.h,v $
                                                                        
// Written: ZHY
// Created: August 2000

// Description: This file contains the class prototype for PressureIndependMultiYield.
//
// What: "@(#) PressureIndependMultiYield.h, revA" 
#ifndef PressureIndependMultiYield_h
#define PressureIndependMultiYield_h

#include <NDMaterial.h>
#include "soil/T2Vector.h"
#include <Matrix.h>

class MultiYieldSurface;

class PressureIndependMultiYield : public NDMaterial
{
public:
     // Initialization constructor
     PressureIndependMultiYield (int tag, 
				 int nd,
				 double rho, 
				 double refShearModul,
				 double refBulkModul,
				 double cohesi,
				 double peakShearStra,
				 double frictionAng = 0.,
				 double refPress = 100, 
				 double pressDependCoe = 0.0,
				 int   numberOfYieldSurf = 20,
				 double * gredu = 0);

     // Default constructor
     PressureIndependMultiYield ();

     // Copy constructor
     PressureIndependMultiYield (const PressureIndependMultiYield &);

     // Destructor: clean up memory storage space.
     virtual ~PressureIndependMultiYield ();

     const char *getClassType(void) const {return "PressureIndependMultiYield";};

     double getRho(void) {return rhox[matN];} ;

     // Sets the values of the trial strain tensor.
     int setTrialStrain (const Vector &strain);

     // Sets the values of the trial strain and strain rate tensors.
     int setTrialStrain(const Vector &v, const Vector &r);

     int setTrialStrainIncr(const Vector &v);
     int setTrialStrainIncr(const Vector &v, const Vector &r);

     // Calculates current tangent stiffness.
     const Matrix &getTangent (void);
     const Matrix &getInitialTangent (void);
        
     void getBackbone (Matrix &);

     // Calculates the corresponding stress increment (rate), for a given strain increment. 
     const Vector &getStress (void);
     const Vector &getStrain (void);
     const Vector &getCommittedStress (void);
     const Vector &getStressToRecord (int numOutput); // Added by Alborz Ghofrani - UW
     const Vector &getCommittedStrain (void);

     // Accepts the current trial strain values as being on the solution path, and updates 
     // all model parameters related to stress/strain states. Return 0 on success.
     int commitState (void);

     // Revert the stress/strain states to the last committed states. Return 0 on success.
     int revertToLastCommit (void);

     int revertToStart(void) {return 0;}

     // Return an exact copy of itself.
     NDMaterial *getCopy (void);

     // Return a copy of itself if "code"="PressureIndependMultiYield", otherwise return null.
     NDMaterial *getCopy (const char *code);

     // Return the string "PressureIndependMultiYield".
     const char *getType (void) const ;

     // Return ndm.
     int getOrder (void) const ;

     int sendSelf(int commitTag, Channel &theChannel);  
     int recvSelf(int commitTag, Channel &theChannel, 
		  FEM_ObjectBroker &theBroker);    

     Response *setResponse (const char **argv, int argc, OPS_Stream &s);
     int getResponse (int responseID, Information &matInformation);
     void Print(OPS_Stream &s, int flag =0);

     //void setCurrentStress(const Vector stress) { currentStress=T2Vector(stress); }
     int setParameter(const char **argv, int argc, Parameter &param);
     int updateParameter(int responseID, Information &eleInformation);	

    // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial stresses.
    friend class PyLiq1;
    friend class TzLiq1;

protected:

private:

	static int matCount;
	static int* loadStagex;  //=0 if elastic; =1 if plastic

  // user supplied
	static int* ndmx;  //num of dimensions (2 or 3)
	static double* rhox;
	static double* frictionAnglex;
	static double* peakShearStrainx;
	static double* refPressurex;
	static double* cohesionx;
	static double* pressDependCoeffx;
	static int*    numOfSurfacesx;

	// internal
	static double* residualPressx;
	static Matrix theTangent;  //classwise member
	int e2p;
	int matN;
	double refShearModulus;
	double refBulkModulus;
	MultiYieldSurface * theSurfaces; // NOTE: surfaces[0] is not used  
	MultiYieldSurface * committedSurfaces;  
	int    activeSurfaceNum;  
	int    committedActiveSurf;
	T2Vector currentStress;
	T2Vector trialStress;
	T2Vector currentStrain;
	T2Vector strainRate;
	static T2Vector subStrainRate;
    double * mGredu;

	void elast2Plast(void);
	// Called by constructor
	void setUpSurfaces(double *);  

	double yieldFunc(const T2Vector & stress, const MultiYieldSurface * surfaces, 
			 int surface_num);

	void deviatorScaling(T2Vector & stress, const MultiYieldSurface * surfaces, 
			     int surfaceNum, int count=0);

	void initSurfaceUpdate(void);

	void paramScaling(void);

	// Return num_strain_subincre
	int setSubStrainRate(void);

	int isLoadReversal(void);

	void getContactStress(T2Vector &contactStress);

	void getSurfaceNormal(const T2Vector & stress, Vector &surfaceNormal); 

	void setTrialStress(T2Vector & stress); 

	double getLoadingFunc(const T2Vector & contact, 
			      const Vector & surfaceNormal,int crossedSurface);

	void stressCorrection(int crossedSurface);
	
	void updateActiveSurface(void);
	
	void updateInnerSurface(void);

	// Return 1 if crossing the active surface; return 0 o/w
	int  isCrossingNextSurface(void);  

};

#endif







