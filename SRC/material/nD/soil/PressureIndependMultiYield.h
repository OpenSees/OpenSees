// $Revision: 1.13 $
// $Date: 2003-02-25 23:33:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PressureIndependMultiYield.h,v $
                                                                        
// Written: ZHY
// Created: August 2000

// Description: This file contains the class prototype for PressureIndependMultiYield.
//
// What: "@(#) PressureIndependMultiYield.h, revA"

#ifndef PressureIndependMultiYield_h
#define PressureIndependMultiYield_h

#include <NDMaterial.h>
#include <MultiYieldSurface.h>
#include <Matrix.h>
#include <Tensor.h>

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

		 double getRho(void) {return rho;} ;
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
     const Vector &getCommittedStrain (void);

     int setTrialStrain (const Tensor &v) {return 0;}
     int setTrialStrain (const Tensor &v, const Tensor &r) {return 0;}
     int setTrialStrainIncr (const Tensor &v) {return 0;}
     int setTrialStrainIncr (const Tensor &v, const Tensor &r) {return 0;}

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

     Response *setResponse (const char **argv, int argc, Information &matInfo);
     int getResponse (int responseID, Information &matInformation);
     void Print(OPS_Stream &s, int flag =0);

     //void setCurrentStress(const Vector stress) { currentStress=T2Vector(stress); }
     int updateParameter(int responseID, Information &eleInformation);	

protected:

private:
	int ndm;  //num of dimensions (2 or 3)
	static int loadStage;  //=0 if elastic; =1 if plastic

  // user supplied
	double rho;
	double refShearModulus;
	double refBulkModulus;
	double frictionAngle;
	double peakShearStrain;
	double refPressure;
	double cohesion;
	double pressDependCoeff;
	int    numOfSurfaces;

	// internal
	double residualPress;
	static Matrix theTangent;  //classwise member
	int e2p;
	MultiYieldSurface * theSurfaces; // NOTE: surfaces[0] is not used  
	MultiYieldSurface * committedSurfaces;  
	int    activeSurfaceNum;  
	int    committedActiveSurf;
	T2Vector currentStress;
	T2Vector trialStress;
	T2Vector currentStrain;
	T2Vector strainRate;
	static T2Vector subStrainRate;

	void elast2Plast(void);
	// Called by constructor
	void setUpSurfaces(double *);  

	double yieldFunc(const T2Vector & stress, const MultiYieldSurface * surfaces, 
			 int surface_num);

	void deviatorScaling(T2Vector & stress, const MultiYieldSurface * surfaces, 
			     int surfaceNum);

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



