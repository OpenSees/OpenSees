// $Revision: 1.1 $
// $Date: 2000-12-19 03:35:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PressureDependMultiYield.h,v $
                                                                        
// Written: ZHY
// Created: August 2000


// Description: This file contains the class prototype for PressureDependMultiYield.
//
// What: "@(#) PressureDependMultiYield.h, revA"

#ifndef PressureDependMultiYield_h
#define PressureDependMultiYield_h

#include <NDMaterial.h>
#include "MultiYieldSurface.h"
#include <Matrix.h>
#include <Tensor.h>

#define MAT_TAG_PressureDependMultiYield		101

class PressureDependMultiYield : public NDMaterial
{
public:
     // Initialization constructor
     PressureDependMultiYield (int tag, 
			              int nd,
                    double refShearModul,
    						  	double refBulkModul,
                    double frictionAng,
										double peakShearStra,
										double refPress,
										double cohesi,
										double pressDependCoe,
                    int   numberOfYieldSurf,
										double phaseTransformAngle, 
                    double contractionParam1,
                    double contractionParam2,
                    double dilationParam1,
                    double dilationParam2,
										double volLimit,
                    double liquefactionParam1,
                    double liquefactionParam2,
                    double liquefactionParam3,
                    double liquefactionParam4,
										double atm);

     // Default constructor
     PressureDependMultiYield ();

     // Copy constructor
     PressureDependMultiYield (const PressureDependMultiYield &);

     // Destructor: clean up memory storage space.
     virtual ~PressureDependMultiYield ();

     // Sets the values of the trial strain tensor.
     int setTrialStrain (const Vector &strain);

		 // Sets the values of the trial strain and strain rate tensors.
		 int setTrialStrain(const Vector &v, const Vector &r);

		 int setTrialStrainIncr(const Vector &v);
     int setTrialStrainIncr(const Vector &v, const Vector &r);

     // Calculates current tangent stiffness.
     const Matrix &getTangent (void);
        
     // Calculates the corresponding stress increment (rate), for a given strain increment. 
     const Vector &getStress (void);
     const Vector &getStrain (void);
	   const Vector &getCommittedStress (void);
     const Vector &getCommittedStrain (void);

		 int setTrialStrain (const Tensor &v) {return 0;}
     int setTrialStrain (const Tensor &v, const Tensor &r) {return 0;}
     int setTrialStrainIncr (const Tensor &v) {return 0;}
     int setTrialStrainIncr (const Tensor &v, const Tensor &r) {return 0;}
     const Tensor &getTangentTensor (void) {Tensor * t=new Tensor; return *t;}
     const Tensor &getStressTensor (void) {Tensor * t=new Tensor; return *t;}
		 const Tensor &getStrainTensor(void) {Tensor * t=new Tensor; return *t;}

     // Accepts the current trial strain values as being on the solution path, and updates 
     // all model parameters related to stress/strain states. Return 0 on success.
     int commitState (void);

     // Revert the stress/strain states to the last committed states. Return 0 on success.
     int revertToLastCommit (void);

		 int revertToStart(void) {return 0;}

     // Return an exact copy of itself.
     NDMaterial *getCopy (void);

     // Return a copy of itself if "code"="PressureDependMultiYield", otherwise return null.
     NDMaterial *getCopy (const char *code);

     // Return the string "PressureDependMultiYield".
     const char *getType (void) const ;

     // Return ndm.
     int getOrder (void) const ;

		 int sendSelf(int commitTag, Channel &theChannel);  
     int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
		 int getResponse (int responseID, Information &matInformation);
     void Print(ostream &s, int flag =0);

		 //void setCurrentStress(const Vector stress) { currentStress=T2Vector(stress); }
     int updateParameter(int responseID, Information &eleInformation);

protected:

private:
  // user supplied 
	int ndm;  //num of dimensions (2 or 3)
	static int loadStage;  //=0 if elastic; =1 or 2 if plastic
	double refShearModulus;
	double refBulkModulus;
	double frictionAngle;
	double peakShearStrain;
	double refPressure;
	double cohesion;
	double pressDependCoeff;
	int    numOfSurfaces;
  double phaseTransfAngle; 
  double contractParam1;
  double contractParam2;
  double dilateParam1;
  double dilateParam2;
	double volumeLimit;
  double liquefyParam1;
  double liquefyParam2;
  double liquefyParam3;
  double liquefyParam4;
	static double AtmoPress;

	// internal
	double residualPress;
	double stressRatioPT;
	static Matrix theTangent;

	int e2p;
  MultiYieldSurface * theSurfaces; // NOTE: surfaces[0] is not used  
  MultiYieldSurface * committedSurfaces;  
  int    activeSurfaceNum;  
	int    committedActiveSurf;
	double modulusFactor;
  T2Vector currentStress;
  T2Vector trialStress;
  T2Vector currentStrain;
  T2Vector strainRate;
  static T2Vector subStrainRate;

	double pressureD;
  T2Vector reversalStress;
	int onPPZ; //=-1 never reach PPZ before; =0 below PPZ; =1 on PPZ; =2 above PPZ
	double strainPTOcta;
	double PPZSize;
	double cumuDilateStrainOcta;
  double maxCumuDilateStrainOcta;
	double cumuTranslateStrainOcta;
	double prePPZStrainOcta;
	double oppoPrePPZStrainOcta;
	static T2Vector trialStrain;
	T2Vector PPZPivot;
	T2Vector PPZCenter;
	T2Vector lockStress;

	double pressureDCommitted;
  T2Vector reversalStressCommitted;
	int onPPZCommitted;
	double PPZSizeCommitted;
	double cumuDilateStrainOctaCommitted;
  double maxCumuDilateStrainOctaCommitted;
	double cumuTranslateStrainOctaCommitted;
	double prePPZStrainOctaCommitted;
	double oppoPrePPZStrainOctaCommitted;
	T2Vector PPZPivotCommitted;
	T2Vector PPZCenterCommitted;
	T2Vector lockStressCommitted;
	static Vector workV;
	static Matrix workM;

	void elast2Plast(void);
	// Called by constructor
  void setUpSurfaces(void);  
  double yieldFunc(const T2Vector & stress, const MultiYieldSurface * surfaces, 
		               int surface_num);
  void deviatorScaling(T2Vector & stress, const MultiYieldSurface * surfaces, 
		                   int surfaceNum);
	void initSurfaceUpdate(void);
  void initStrainUpdate(void);

  // Return num_strain_subincre
  int setSubStrainRate(void);
	int isLoadReversal(void);
	T2Vector getContactStress(void);
  T2Vector getSurfaceNormal(const T2Vector & stress); 
	double getModulusFactor(T2Vector & stress);
	void updatePPZ(const T2Vector & stress);
  void PPZTranslation(const T2Vector & contactStress);
	double getPPZLimits(int which, const T2Vector & contactStress);
	double getPlasticPotential(const T2Vector & stress, const T2Vector & surfaceNormal);
  void setTrialStress(T2Vector & stress); 
  double getLoadingFunc(const T2Vector & contact, const T2Vector & surfaceNormal,
		                    double plasticPotential,int crossedSurface);
	//return 1 if stress locked; o/w return 0.
  int stressCorrection(int crossedSurface);
  void updateActiveSurface(void);
  void updateInnerSurface(void);

  // Return 1 if crossing the active surface; return 0 o/w
  int  isCrossingNextSurface(void);  

};

#endif
