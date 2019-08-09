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

//Modified Ibarra-Medina-Krawinkler with Bilin Hysteretic Response

//**********************************************************************                                                                     
// Code Developed by: Ahmed Elkady
// Postdoctoral Researcher, EPFL, Switzerland
// Last Updated: August 29th 2012
//**********************************************************************

#ifndef IMKBilin_h
#define IMKBilin_h

#include <UniaxialMaterial.h>

class IMKBilin : public UniaxialMaterial
{
public:
	IMKBilin(int tag, double Ke,
		double Theta_p_pos0, double Theta_pc_pos0, double Theta_u_pos0, double Mpe_pos0, double MmaxMpe_pos0, double ResM_pos0,
		double Theta_p_neg0, double Theta_pc_neg0, double Theta_u_neg0, double Mpe_neg0, double MmaxMpe_neg0, double ResM_neg0,
		double LAMBDA_S, double LAMBDA_C, double LAMBDA_K, double c_S, double c_C, double c_K, double D_pos, double D_neg);
	IMKBilin();
	~IMKBilin();
	const char *getClassType(void) const { return "IMKBilin"; };
	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);
	UniaxialMaterial *getCopy(void);
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	void Print(OPS_Stream &s, int flag = 0);


/* 	double r8_max(double x, double y);
	double r8_min(double x, double y);
	int i4_max(int i1, int i2);
	double pchst(double arg1, double arg2);
	int chfev(double x1, double x2, double f1, double f2, double d1, double d2, int ne, double xe[], double fe[], int next[]);
	void spline_pchip_set(int n, double x[], double f[], double d[]);
	void spline_pchip_val(int n, double x[], double f[], double d[], int ne, double xe[], double fe[]);
 */
protected:

private:
	//my functions


	//Fixed input material parameters 
	double Ke;
	double Theta_p_pos0;
	double Theta_pc_pos0;
	double Theta_u_pos0;
	double Mpe_pos0;
	double MmaxMpe_pos0;
	double ResM_pos0;
	double Theta_p_neg0;
	double Theta_pc_neg0;
	double Theta_u_neg0;
	double Mpe_neg0;
	double MmaxMpe_neg0;
	double ResM_neg0;
	double LAMBDA_S;
	double LAMBDA_C;
	double LAMBDA_K;
	double c_S;
	double c_C;
	double c_K;
	double D_pos;
	double D_neg;
	double n;
	double Roffset;
	double LAMBDA_F;
	double c_F;

	//State variables 
	double U, cU;

	//History variables 
	double 			Mr_pos0, cMr_pos0;
	double 			Mr_neg0, cMr_neg0;
	double 	 Theta_y_pos0;
	double    Theta_max_pos0;
	double 	 slope_p_pos0;
	double 	slope_pc_pos0;
	double 		Mmax_pos0;
	double   MpeProject_pos0;
	double  MmaxProject_pos0;

	double 	 Theta_y_neg0;
	double    Theta_max_neg0;
	double 	 slope_p_neg0;
	double 	slope_pc_neg0;
	double 		Mmax_neg0;
	double   MpeProject_neg0;
	double  MmaxProject_neg0;

	double     Ref_Energy_S;
	double     Ref_Energy_C;
	double     Ref_Energy_K;
	double     Ref_Energy_F;

	double    			 K_j_1, cK_j_1;
	double     Theta_y_pos_j_1, cTheta_y_pos_j_1;
	double   Theta_max_pos_j_1, cTheta_max_pos_j_1;
	double     slope_p_pos_j_1, cslope_p_pos_j_1;
	double    slope_pc_pos_j_1, cslope_pc_pos_j_1;
	double         Mpe_pos_j_1, cMpe_pos_j_1;
	double  MpeProject_pos_j_1, cMpeProject_pos_j_1;
	double        Mmax_pos_j_1, cMmax_pos_j_1;
	double MmaxProject_pos_j_1, cMmaxProject_pos_j_1;

	double     Theta_y_neg_j_1, cTheta_y_neg_j_1;
	double   Theta_max_neg_j_1, cTheta_max_neg_j_1;
	double     slope_p_neg_j_1, cslope_p_neg_j_1;
	double    slope_pc_neg_j_1, cslope_pc_neg_j_1;
	double         Mpe_neg_j_1, cMpe_neg_j_1;
	double  MpeProject_neg_j_1, cMpeProject_neg_j_1;
	double        Mmax_neg_j_1, cMmax_neg_j_1;
	double MmaxProject_neg_j_1, cMmaxProject_neg_j_1;

	double Ri, cRi;
	double Mi, cMi;
	double Di, cDi;
	double Ri_1, cRi_1;
	double Mi_1, cMi_1;
	double Di_1, cDi_1;

	double   beta_S_j_1, cbeta_S_j_1;
	double   beta_C_j_1, cbeta_C_j_1;
	double   beta_K_j_1, cbeta_K_j_1;
	double 	 beta_F_j_1, cbeta_F_j_1;

	double Excursion_Flag, cExcursion_Flag;
	double  Reversal_Flag, cReversal_Flag;
	double 	   Yield_Flag, cYield_Flag;
	double   Fail_FlagPos, cFail_FlagPos;
	double   Fail_FlagNeg, cFail_FlagNeg;
	double 	   Mrpos_Flag, cMrpos_Flag;
	double 	   Mrneg_Flag, cMrneg_Flag;
	double 	  Energy_Flag, cEnergy_Flag;

	double Energy_Excrsni_1, cEnergy_Excrsni_1;
	double 	  Energy_Excrsn, cEnergy_Excrsn;
	double  	 Energy_Rev, cEnergy_Rev;
	double     Energy_total, cEnergy_total;

	double Rreversal, cRreversal;
	double Mreversal, cMreversal;
	double TangentK, cTangentK;

	double beta_Sx, cbeta_Sx;
	double slope_p_x, cslope_p_x;
	double Mpe_x, cMpe_x;
	double Theta_y_x, cTheta_y_x;
	double MpeProject_x, cMpeProject_x;

};


#endif

//////////////////////////// Variables Definitions /////////////////////////////
/*
Ke 					Initial elastic stiffness
Theta_p_pos0 		Initial pre-capping plastic rotation in the +ve loading direction
Theta_pc_pos0   	Initial post-capping plastic rotation in the +ve loading direction
Theta_u_pos0 		Ultimate rotation in the +ve loading direction
Mpe_pos0 			Initial effective plastic moment in the +ve loading direction
MmaxMpe_pos0 		Initial maximum-to-effective plastic moment ratio in the +ve loading direction
ResM_pos0 			Residual moment in the +ve loading direction
Theta_p_neg0 		Initial pre-capping plastic rotation in the -ve loading direction
Theta_pc_neg0   	Initial post-capping plastic rotation in the -ve loading direction
Theta_u_neg0    	Ultimate rotation in the -ve loading direction
Mpe_neg0        	Initial effective plastic moment in the -ve loading direction
MmaxMpe_neg0    	Initial maximum-to-effective plastic moment ratio in the -ve loading direction
ResM_neg0       	Residual moment in the -ve loading direction
LAMBDA_S 			Cyclic deterioration parameter for strength deterioration 
LAMBDA_C 			Cyclic deterioration parameter for post-capping strength deterioration
LAMBDA_K 			Cyclic deterioration parameter for unloading stiffness deterioration 
c_S 				Rate of strength deterioration.
c_C 				Rate of post-capping strength deterioration.
c_K 				Rate of unloading stiffness deterioration.
D_pos 				Rate of cyclic deterioration in the +ve loading direction
D_neg 				Rate of cyclic deterioration in the -ve loading direction
n					Paramter identifying the offset rotation on the unloading side
Roffset				Offset rotation identifying the Smooth Transition region
LAMBDA_F			Cyclic deterioration parameter for Smooth Transition deterioration 
c_F					Cyclic deterioration parameter for Smooth Transition deterioration
Ri 					Rotation at current step
Mi 					Moment at current step 
Di 					Rotation Direction at current step 
Ri_1  				Rotation at previous step
Mi_1            	Moment at previous step 
Di_1            	Rotation Direction at previous step
Rreversal 			Rotation at direction reversal points
Mreversal 			Moment at direction reversal points
TangentK 			Tangent stiffness
K_j_1 				Unloading stiffness in the previous excursion
Theta_y_pos_j_1 	Yielding rotation in the previous +ve excursion
Theta_max_pos_j_1 	Capping point rotation in the previous +ve excursion
slope_p_pos_j_1 	Pre-capping slope in the previous +ve excursion
slope_pc_pos_j_1 	Post-capping slope in the previous +ve excursion
Mpe_pos_j_1 		Effective plastic moment in the previous +ve excursion
MpeProject_pos_j_1  Projected effective plastic moment in the previous +ve excursion
Mmax_pos_j_1        Maximum moment in the previous +ve excursion
MmaxProject_pos_j_1 Projected maximum  moment in the previous +ve excursion
Theta_y_neg_j_1     Yielding rotation in the previous -ve excursion
Theta_max_neg_j_1   Capping point rotation in the previous -ve excursion
slope_p_neg_j_1     Pre-capping slope in the previous -ve excursion
slope_pc_neg_j_1    Post-capping slope in the previous -ve excursion
Mpe_neg_j_1         Effective plastic moment in the previous -ve excursion
MpeProject_neg_j_1  Projected effective plastic moment in the previous -ve 
Mmax_neg_j_1        Maximum moment in the previous -ve excursion
MmaxProject_neg_j_1 Projected maximum  moment in the previous -ve excursion
beta_S_j_1    
beta_C_j_1    
beta_K_j_1    
beta_F_j_1 	  
Ref_Energy_S    	Refernence energy for strength deterioration
Ref_Energy_C    	Refernence energy for post-capping strength deterioration
Ref_Energy_K    	Refernence energy for unloading stiffness deterioration
Ref_Energy_F    	Refernence energy for Smooth Transition deterioration
Excursion_Flag 		Flag for Excursion occurance (i.e., crossing the x-axis)
Reversal_Flag 		Flag for Loading direction reversal occurance
Yield_Flag 			Flag for Yielding occurance
Fail_FlagPos 		Flag for reaching the ultimate rotation in the +ve loading direction
Fail_FlagNeg 		Flag for reaching the ultimate rotation in the -ve loading direction
Mrpos_Flag 			Flag for reaching the residual moment in the +ve loading direction
Mrneg_Flag 			Flag for reaching the residual moment in the -ve loading direction
Energy_Flag 		Flag for reaching the reference energy
Energy_Excrsni_1	Dissipated energy in previous excursion
Energy_Excrsn 		Dissipated energy in current excursion
Energy_Rev 			Total dissipated energy till previous load reversal point
Energy_total 		Total dissipated energy till current step
*/