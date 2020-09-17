// $Revision: 1.1 $
// $Date: 2009-07-23 23:44:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/MultiYieldSurfaceClay.h,v $
                                                                        
// Consistent tangent and Sensitivity Written: Quan Gu UCSD
// Created: Jul 2009
//
// refer to "Finite element response sensitivity analysis of multi-yield-surface J2 
// plasticity model by direct differentiation method", Quan Gu, Joel P. Conte, Ahmed Elgamal
// and Zhaohui Yangc, Computer Methods in Applied Mechanics and Engineering
// Volume 198, Issues 30-32, 1 June 2009, Pages 2272-2285 

// Description: This file contains the class prototype for MultiYieldSurfaceClay.
//
// What: "@(#) MultiYieldSurfaceClay.h, revA"

#ifndef MultiYieldSurfaceClay_h
#define MultiYieldSurfaceClay_h

#include <NDMaterial.h>
#include <Matrix.h>
#include "soil/T2Vector.h"
class MultiYieldSurface;

#define ND_TAG_MultiYieldSurfaceClay   10284765

class MultiYieldSurfaceClay : public NDMaterial
{
public:
     // Initialization constructor
     MultiYieldSurfaceClay (int tag, 
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
     MultiYieldSurfaceClay ();

     // Copy constructor
     MultiYieldSurfaceClay (const MultiYieldSurfaceClay &);

     // Destructor: clean up memory storage space.
     virtual ~MultiYieldSurfaceClay ();

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
     const Vector &getCommittedStrain (void);

     // Accepts the current trial strain values as being on the solution path, and updates 
     // all model parameters related to stress/strain states. Return 0 on success.
     int commitState (void);

     // Revert the stress/strain states to the last committed states. Return 0 on success.
     int revertToLastCommit (void);

     int revertToStart(void) ;

     // Return an exact copy of itself.
     NDMaterial *getCopy (void);

     // Return a copy of itself if "code"="MultiYieldSurfaceClay", otherwise return null.
     NDMaterial *getCopy (const char *code);

     // Return the string "MultiYieldSurfaceClay".
     const char *getType (void) const ;

     // Return ndm.
     int getOrder (void) const ;

     int sendSelf(int commitTag, Channel &theChannel);  
     int recvSelf(int commitTag, Channel &theChannel, 
		  FEM_ObjectBroker &theBroker);    

     Response *setResponse (const char **argv, int argc, OPS_Stream &theOutput);
     int getResponse (int responseID, Information &matInformation);
     void Print(OPS_Stream &s, int flag =0);

     //void setCurrentStress(const Vector stress) { currentStress=T2Vector(stress); }
//	 int updateParameter(int responseID, Information &eleInformation,int Yang);
	
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

	////////////////////add sensitivity ////////////////////////
public:
    int            setParameter             (const char **argv, int argc, Parameter &param);
//    int            setParameter             (const char **argv, int argc, Information &info);
    int            updateParameter          (int parameterID, Information &info);
	int            activateParameter        (int parameterID);
	const Vector & getStressSensitivity     (int gradNumber, bool conditional);
	int            commitSensitivity        (const Vector & strainGradient, int gradNumber, int numGrads);


    void       setUpSurfacesSensitivity(int GradNumber);
	void       setTrialStressSensitivity(T2Vector & stress,T2Vector & dStress);
    void       updateInnerSurfaceSensitivity(void);
    int        setSubStrainRateSensitivity(void);
	void       getContactStressSensitivity(T2Vector &contactStress);

    void       getSurfaceNormalSensitivity(const T2Vector & stress,const T2Vector & dStress1, 
	                                     Vector &surfaceNormal,Vector &dSurfaceNormal);
	double     getLoadingFuncSensitivity(const T2Vector & contactStress, 
									 const Vector & surfaceNormal,const Vector & dSurfaceNormal,
									 int crossedSurface);
	void       updateActiveSurfaceSensitivity(void);
	void       stressCorrectionSensitivity(int crossedSurface);
	int        isSurfacesSensitivitySetUp(int gradNumber);
    void       setSurfacesSensitivitySetUpMark(int gradNumber)	;

	    // new methods for recorder requested by Quan .. Feb 1 05
	const  Vector &getCommittedStressSensitivity(int GradientNumber);
	const  Vector &getCommittedStrainSensitivity(int GradientNumber); 
private:	
	// Sensitivity history variables
	Matrix *SHVs;
    // Parameter identification
	int parameterID;
    // define more sensitivity to store data
  	int myNumGrads;
	int gradNumber;
	double dVolume;
	double dLoadingFunc;
	int debugMarks;               // NONclasswide int for debug only

	static T2Vector dCurrentStress;
	static T2Vector dTrialStress;
	static T2Vector dCurrentStrain;
	static T2Vector dSubStrainRate;
	static T2Vector dStrainRate;
	static T2Vector dContactStress;

// uncommitted conditional sensitivity
	double * dMultiSurfaceCenter;

// committed unconditional sensitivity

	double * dCommittedMultiSurfaceSize;
//	double * dCommittedMultiSurfaceElastPlastModul;
	double * dCommittedMultiSurfacePlastModul;
	double * dCommittedMultiSurfaceCenter;
	int    * surfacesSensitivityMark;
	

////////////////////// end sensitivity ///////////////////////


//------------- for consistent tangent ----------------
private:
	static Matrix dTrialStressdStrain;   //classwide matrix
	Matrix consistentTangent;
	static Matrix dContactStressdStrain;  // classwide matrix
	static Matrix dSurfaceNormaldStrain;  // classwide matrix
	static Vector dXdStrain;              // classwide Vector

	static Vector     temp6;          // classwide Vector
	static Vector     temp;           // classwide Vector
	static Vector     devia;          // classwide Vector



};

#endif



