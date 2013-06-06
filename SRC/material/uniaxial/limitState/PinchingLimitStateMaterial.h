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
                                                                  
#ifndef PinchingLimitStateMaterial_h
#define PinchingLimitStateMaterial_h

// Written: MRL 
//
// Description: This file contains the class implementation for PinchingLimitStateMaterial.
//
// What: "@(#) PinchingLimitStateMaterial.h, revN/C"

#include <UniaxialMaterial.h>
#include <LimitCurve.h>
#include <Matrix.h>

class Element;
class Domain;
class DomainComponent;
class Node;


class PinchingLimitStateMaterial : public UniaxialMaterial
{
	public:
		PinchingLimitStateMaterial(int matTag, 
		int nodeT, int nodeB, int drftAx, double Kelas, int crvTyp, int crvTag,
		double YpinchUPN, double YpinchRPN,    double XpinchRPN,
		double YpinchUNP, double YpinchRNP,    double XpinchRNP,
		double dmgStrsLimE, double dmgDispMax,
		double dmgE1, double dmgE2, double dmgE3, double dmgE4, double dmgELim,
		double dmgU1, double dmgU2, double dmgU3, double dmgU4, double dmgULim,
		double dmgR1, double dmgR2, double dmgR3, double dmgR4,	double dmgRLim, double dmgRCyc,
		double dmgS1, double dmgS2, double dmgS3, double dmgS4, double dmgSLim, double dmgSCyc, 
		int eleTag,	double b, double d, double h, double a, double st, double As, double Acc,
		double ld, double db, double rhot, double fc,double fy, double fyt, 
		Domain *theDomain, Node *theNodeT, Node *theNodeB, LimitCurve &theCurve, Element *theElement);			

		PinchingLimitStateMaterial();    

		~PinchingLimitStateMaterial();
		
		int setTrialStrain(double strain, double strainRate = 0.0); 
		double getStrain(void);          
		double getStress(void);
		double getTangent(void);
		double getInitialTangent(void) {return E1;};

		int commitState(void);
		int revertToLastCommit(void);  
		int revertToStart(void);	

		UniaxialMaterial *getCopy(void);
    
		int sendSelf(int commitTag, Channel &theChannel);  
		int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
		void Print(OPS_Stream &s, int flag =0);

	protected:
    
	private:

		//Added Member Functions
		void updateDamageE(void);	
		void updateDamageS(void);	
		void updateDamageR(void);	
		void updateEnergy(void);					
		void defineBackbone(void);	
		int getStateFlag(void);			
		void definePinchingPN(void);	
		void definePinchingNP(void);
		double getFlexDisp(void);	
		double getFlexShift(void);
		void checkEnvelope(void);	
		void defineTargetVars(void);
		double getAxialForce(void);	
		void defineE1(void);

		// Input Variables
		int nodeTop;	// Top node of the column flexural section
		int nodeBot;	// Bottom node of the column flexural section
		int driftAxis;	// Axis of flexural deformation
						// driftAxis = 1 -- Drift along the x-axis
						// driftAxis = 2 -- Drift along the y-axis
						// driftAxis = 3 -- Drift along the z-axis
		double E1;		// Initial elastic stiffness
		int curveType;	// type of limit curve = 0 for no curve
						//						 1 for axial curve
						//						 2 for shear curve
		int curveTag;	// unique tag identifying associated limit curve
		int eleTag;		// tag for associated beam column element
		double YpinchUnloadPN;						// Factor of load at reversal to set the stress component of the unloading breakpoint in the Pos->Neg direction (strain breakpoint set by unloading slope)
		double XpinchReloadPN, YpinchReloadPN;		// Factor of load at reversal to set the strain and stress components of the reloading breakpoint in the Pos->Neg direction
		double YpinchUnloadNP;						// Factor of load at reversal to set the stress component of the unloading breakpoint in the Neg->Pos direction (strain breakpoint set by unloading slope)
		double XpinchReloadNP, YpinchReloadNP;		// Factor of load at reversal to set the strain and stress components of the reloading breakpoint in the Neg->Pos direction
		double dmgStressLimE;		// stress at proportional limit
		double dmgDeflMax;			// Maximum estimated displacement of limit spring at column collapse should compensate for member deflection if using tip deflection
		double dmgElastic1, dmgElastic2, dmgElastic3, dmgElastic4, dmgElasticLim;						// Damage parameters for elastic stiffness degredation
		double dmgUnload1, dmgUnload2, dmgUnload3, dmgUnload4, dmgUnloadLim;							// Damage parameters for unloading stiffness degredation
		double dmgReload1, dmgReload2, dmgReload3, dmgReload4, dmgReloadLim, dmgReloadCyclic;			// Damage parameters for reloading stiffness degredation, dmgReloadCyclic is a factor of previous load at unloading
		double dmgStrength1, dmgStrength2, dmgStrength3, dmgStrength4, dmgStrengthLim, dmgStrengthCyclic;// Damage parameters for strength degredation

		Domain *theDomain;		
		Node *theNodeT;			
		Node *theNodeB;		
		LimitCurve *theCurve;	
		Element *theElement;

		double b;			// column width (in)
		double d;			// column effective depth (in)
		double h;			// column cross section height (in)
		double a;			// column shear span length (in)
		double st;			// transverse reinforcement spacing (in)
		double As;			// area of longitudinal steel bars in column section (in^2)
		double Acc;			// gross confined concrete area bounded by the transverse reinforcement in column section (in^2)
		double ld;			// development length of longitudinal bars using ACI 318-08 Equations 12-1 & 12-2 (in)
		double db;			// diameter of longitudinial bars (in)
		double rhot;		// transverse reinforcement ratio Ast/(st*b)
		double fc;			// concrete strength (ksi)
		double fy;			// longitudinal steel yield strength (ksi)
		double fyt;			// transverse steel yield strength (ksi)


		// Trial Variables
		double Tstress;
		double Tstrain;
		double Ttangent;
		double TstrainRate;	
		double TstrainMax;		// Maximum strain
		double TstrainMin;		// Minimum strain
		double Tdu;				// Change in strain
		double TpinchStressUnloadPN, TpinchStrainUnloadPN;		// Pinching stress and strain breakpoint during unloading from positive to negative direction
		double TpinchStressReloadPN, TpinchStrainReloadPN;		// Pinching stress and strain stress breakpoint during reloading from positive to negative direction
		double TpinchStressUnloadNP, TpinchStrainUnloadNP;		// Pinching stress and strain stress breakpoint during unloading from negative to positive direction
		double TpinchStressReloadNP, TpinchStrainReloadNP;		// Pinching stress and strain stress breakpoint during reloading from negative to positive direction
		double TbUnloadPN;				// y-intercept of unloading curve Pos->Neg
		double TbUnloadNP;				// y-intercept of unloading curve Neg->Pos
		int TstateFlag;					// Flag indicating if the limit state surface has been reached
		double TdmgElasticE;			// Elastic damage stiffness
		double TdmgReloadE;				// Reloading damage stiffness
		double TbKdegDmg;				// Damaged y intercept of degrading slope after failure for strength reduction
		double TpinchSlopePN;			// Slope of pinching curve from Pos->Neg
		double TpinchInterceptPN;		// y-intercept of pinching curve Pos->Neg
		double TreloadInterceptPN;		// y-intercept of reloading curve Pos->Neg
		double TpinchSlopeNP;			// Slope of pinching curve from Neg->Pos
		double TpinchInterceptNP;		// y-intercept of pinching curve Neg->Pos
		double TreloadInterceptNP;		// y-intercept of reloading curve Neg->Pos
		double Tenergy;					// current area calculated during loading on envelope curve
		double TstrainFresKdegDmg;		// strain at the point where degrading slope intersects the residual force after failure updated for damage
		double TstrainShearFailureDmg; 	// strain at the point where shear failure occurs updated for damage
		double TbReloadAfterUnloadPN;	// y-intercept of reloading curve Pos->Neg when negative envelope curve is not reached during reloading in the Neg->Pos direction
		double TbReloadAfterUnloadNP;	// y-intercept of reloading curve Neg->Pos when positive envelope curve is not reached during reloading in the Pos->Neg direction
		double TstrainFlex;				// the current flexural displacement from the last iteration
		double TstrainGlobal;			// global displacement
		double TenergyE;				// elastic energy

		// Converged Material History parameters
		double Cstress;
		double Cstrain;
		double Ctangent;	
		double CstrainRate;
		double CstrainMax;		// Maximum strain
		double CstrainMin;		// Minimum strain
		double Cdu;				// Change in strain
		double CpinchStressUnloadPN, CpinchStrainUnloadPN;		// Pinching stress and strain breakpoint during unloading from positive to negative direction
		double CpinchStressReloadPN, CpinchStrainReloadPN;		// Pinching stress and strain stress breakpoint during reloading from positive to negative direction
		double CpinchStressUnloadNP, CpinchStrainUnloadNP;		// Pinching stress and strain stress breakpoint during unloading from negative to positive direction
		double CpinchStressReloadNP, CpinchStrainReloadNP;		// Pinching stress and strain stress breakpoint during reloading from negative to positive direction
		double CbUnloadPN;		// y-intercept of unloading curve Pos->Neg
		double CbUnloadNP;		// y-intercept of unloading curve Neg->Pos
		int CstateFlag;			// Flag indicating if the limit state surface has been reached
								// (not reached = 0, reached for first time on positive envelope = 1, 
								// reached for first time on negative envelope = 2, on positive envelope = 3,
								// on negative envelope = 4, reloading pos->neg = 5,  reloading pos->neg = 5)
		double CdmgElasticE;	// Elastic Damage Stiffness
		double CdmgReloadE;		// Reloading damage stiffness
		double CbKdegDmg;		// Damaged y intercept of degrading slope after failure for strength reduction
		double CpinchSlopePN;	// Slope of pinching curve from Pos->Neg
		double CpinchSlopeNP;	// Slope of pinching curve from Neg->Pos
		double CpinchInterceptPN;		// y-intercept of pinching curve Pos->Neg
		double CpinchInterceptNP;		// y-intercept of pinching curve Neg->Pos
		double CreloadInterceptPN;		// y-intercept of reloading curve Pos->Neg
		double CreloadInterceptNP;		// y-intercept of reloading curve Neg->Pos
		double Cenergy;					// current area calculated during loading on envelope curve
		double CstrainFresKdegDmg;		// strain at the point where degrading slope intersects the residual force after failure updated for damage
		double CstrainShearFailureDmg; 	// strain at the point where shear failure occurs updated for damage
		double CbReloadAfterUnloadPN;	// y-intercept of reloading curve Pos->Neg when negative envelope curve is not reached during reloading in the Neg->Pos direction
		double CbReloadAfterUnloadNP;	// y-intercept of reloading curve Neg->Pos when positive envelope curve is not reached during reloading in the Pos->Neg direction
		double CstrainFlex;				// the current flexural displacement from the last iteration
		double CstrainGlobal;			// global displacement
		double CenergyE;				// elastic energy

		// Static Variables
		double Kdeg;					// degrading slope after failure
		double bKdeg;					// y intercept of degrading slope after failure
		double Fres;					// residual force after failure
		double strainShearFailure;		// strain at the point where shear failure occurs
		double stressShearFailure;		// stress at the point where shear failure occurs
		double InelastMonoEnergy;		// Inelastic monotonic energy under envelope curve for damage computations
		double strainFlexRev;			// flexural strain at loading reversal
		double strainUn;				// strain at unloading from the envelope
		double stressUn;				// stress at unloading from the envelope
		double strainReload;			// strain of material when reloading is indicated
		double stressReload;			// stress of material when reloading is indicated
		double Ereload;					// the initial slope during reloading from the current point to the point on the envelope where unloading initiated
		double slopeFlexPred;			// predicted flexural reloading slope to match displacement at the load reversal point
		double interceptFlexPred;		// predicted flexural intercept of predicted flexural reloading slope
		int resFlag;					// flag to indicate when the residual strength capacity has been reached 0 before residual strength, 1 after residual strength
		int strainFlexRevFlag;
		double stressUnDmg;
		double strainUnDmg;
		double slopeGlobalEnv;
		double interceptGlobalEnv;
		double strainGlobalFresKdeg;
		double strainFlexRevDmg;
		int countGlobalEnv;	//used to allow the stiffness to stabalize before calculating envelope
};

#endif



