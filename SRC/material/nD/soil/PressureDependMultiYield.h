// $Revision: 1.22 $
// $Date: 2007-03-30 01:52:10 $
//>>>>>>> 1.7

// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PressureDependMultiYield.h,v $
// Written: ZHY
// Created: August 2000


// Description: This file contains the class prototype for PressureDependMultiYield.
//
// What: "@(#) PressureDependMultiYield.h, revA"

#ifndef PressureDependMultiYield_h
#define PressureDependMultiYield_h

#include <NDMaterial.h>
#include "soil/T2Vector.h"
#include <Matrix.h>

class MultiYieldSurface;

class PressureDependMultiYield : public NDMaterial
{
public:
     // Initialization constructor
     PressureDependMultiYield (int tag, 
			       int nd,
						 double rho,
			       double refShearModul,
			       double refBulkModul,
			       double frictionAng,
			       double peakShearStra,
			       double refPress,
			       double pressDependCoe,
			       double phaseTransformAngle, 
			       double contractionParam1,
			       double dilationParam1,
			       double dilationParam2,
			       double liquefactionParam1,
			       double liquefactionParam2,
			       double liquefactionParam4,
			       int   numberOfYieldSurf = 20,
				   double * gredu = 0,
		           double e = 0.6,
			       double volLimit1 = 0.9,
			       double volLimit2 = 0.02,
			       double volLimit3 = 0.7,
			       double atm = 101.,
				   double cohesi = 0.1,
				   double hv = 0.,
				   double pv = 1.);

     // Default constructor
     PressureDependMultiYield ();

     // Copy constructor
     PressureDependMultiYield (const PressureDependMultiYield &);

     // Destructor: clean up memory storage space.
     virtual ~PressureDependMultiYield ();

     const char *getClassType(void) const {return "PressureDependMultiYield";};     
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

     // Return a copy of itself if "code"="PressureDependMultiYield", otherwise return null.
     NDMaterial *getCopy (const char *code);

     // Return the string "PressureDependMultiYield".
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
    // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial stresses.    friend class PyLiq1;    friend class TzLiq1;
protected:

private:
  // user supplied 
	 static int matCount;
     static int* ndmx;  //num of dimensions (2 or 3)
     static int* loadStagex;  //=0 if elastic; =1 or 2 if plastic
     static double* rhox;  //mass density
     static double* refShearModulusx;
     static double* refBulkModulusx;
     static double* frictionAnglex;
     static double* peakShearStrainx;
     static double* refPressurex;
     static double* cohesionx;
     static double* pressDependCoeffx;
     static int*    numOfSurfacesx;
     static double* phaseTransfAnglex; 
     static double* contractParam1x;
     static double* dilateParam1x;
     static double* dilateParam2x;
     static double* liquefyParam1x;
     static double* liquefyParam2x;
     static double* liquefyParam4x;
     static double* einitx;    //initial void ratio
     static double* volLimit1x;
     static double* volLimit2x;
     static double* volLimit3x;
     static double pAtm;
	 static double* Hvx;
	 static double* Pvx;

     // internal
     static double* residualPressx;
     static double* stressRatioPTx;
     static Matrix theTangent;
     
	 int matN;
     int e2p;
     MultiYieldSurface * theSurfaces; // NOTE: surfaces[0] is not used  
     MultiYieldSurface * committedSurfaces;  
     int    activeSurfaceNum;  
     int    committedActiveSurf;
     double modulusFactor;
	 double initPress;
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
     static Vector workV6;
     static T2Vector workT2V;
	 double maxPress;
     
     void elast2Plast(void);
     // Called by constructor
     void setUpSurfaces(double *);  
     double yieldFunc(const T2Vector & stress, const MultiYieldSurface * surfaces, 
		      int surface_num);
     void deviatorScaling(T2Vector & stress, const MultiYieldSurface * surfaces, 
			  int surfaceNum);
     void initSurfaceUpdate(void);
     void initStrainUpdate(void);

     // Return num_strain_subincre
     int setSubStrainRate(void);
     int isLoadReversal(const T2Vector &);
     int isCriticalState(const T2Vector & stress);
     void getContactStress(T2Vector &contactStress);
     void getSurfaceNormal(const T2Vector & stress, T2Vector &normal); 
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
