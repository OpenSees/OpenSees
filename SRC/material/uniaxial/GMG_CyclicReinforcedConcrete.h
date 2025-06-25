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

// Written: Rasool Ghorbani, University of Texas at San Antonio
// Created: February 2019
// Description: This file contains the class implementation for GMG_CyclicReinforcedConcrete material Model
// Refer to: <Ghorbani, R., A. Suselo, S. Gendy, A. Matamoros and W. Ghannoum (2022). "Uniaxial model for simulating the cyclic behavior of reinforced concrete members." Earthquake Engineering & Structural Dynamics.>
// Refer to: <Ghorbani, R. (2022). Computational Framework for Decision-Oriented Reinforced Concrete Column Simulation Capabilities. Ph.D., The University of Texas at San Antonio.>

#ifndef GMG_CyclicReinforcedConcrete_h
#define GMG_CyclicReinforcedConcrete_h

#include <UniaxialMaterial.h>

class GMG_CyclicReinforcedConcrete : public UniaxialMaterial
{
public:
	GMG_CyclicReinforcedConcrete(int tag, double Kepos, double Keneg,
		double fypos, double fyneg, double fcappos, double fcapneg,
		double dcappos, double dcapneg, 
		double Kdegpos, double Kdegneg,
		double frespos, double fresneg, 
		double delUpos, double delUneg,
		double alpha_Er_Post_Yielding, double beta_Er_Post_Yielding,
		double alpha_Er_Post_Capping, double beta_Er_Post_Capping,
		double Er_Post_Yielding, double Er_Post_Capping,        
		double Kun_Post_Yielding, double Kun_Post_Capping,  
		double Kr_Post_Yielding, double Kr_Post_Capping,     
		double delta_ratio_max_hard, double Ref_Energy_Coe,
		double C1, double C2, double C3, 
		double solpe_post_yielding, double solpe_post_capping);

	GMG_CyclicReinforcedConcrete();
	~GMG_CyclicReinforcedConcrete();

	const char *getClassType(void) const { return "GMG_CyclicReinforcedConcrete"; };

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStrainRate(void);
	double getStress(void);

	double getTangent(void);
	double getInitialTangent(void);
	double getDampTangent(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	UniaxialMaterial *getCopy(void);

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);

	void Print(OPS_Stream &s, int flag = 0);


protected:

private:

	void splineparam(double MtoRref, double dpeakpos, double dcappose, double dpeakneg, double dcapneg);

	void spline_curve(double dirtag, double P0x, double P0y, double P4x, double P4y, double Krel, double Kun, double E, double x);

	//--- MONOTONIC Backbone
	double Kepos, Keneg, fypos, fyneg, fcappos, fcapneg, dcappos, dcapneg;
	double Kdegpos, Kdegneg, frespos, fresneg, delUpos, delUneg;

	//Spline-Curve properties 
	double alpha_Er_Post_Yielding, beta_Er_Post_Yielding;
	double alpha_Er_Post_Capping, beta_Er_Post_Capping;
	double Er_Post_Yielding, Er_Post_Capping;
	double Kun_Post_Yielding, Kun_Post_Capping;
	double Kr_Post_Yielding, Kr_Post_Capping;

	//--- DAMAGE Properties
	double C1, C2, C3, solpe_post_yielding, solpe_post_capping;
	double delta_ratio_max_hard, Ref_Energy_Coe;
	//double solpe_post_yielding, Csolpe_post_yielding;

	// State Variables
	double d;     // Deformation
	double f;     // Force
	double ek;    // Tangent Stiffness or slope
	double esel;  // Secant to elastic slope
	double MtoRref, KrR, Krel, KuR, Kun, ER, E;
	double x;
	int dirtag;
	double fspl, fderspl, dderspl;
	double din, fin;
	double fouterPN, fouterNP;
	double ekouterPN, ekouterNP;
	double finner, ekinner;
	double s, cwct, cwcc, Epos, Eneg, CEpos, CEneg;

	// Trial and Committeed State Variables
	double Cstrain;         // Committed Strain
	double Cstress;         // Committed Stress
	double Cek;				// Committed Tangent
	double InCycFac, CInCycFac; // In-cycle Factor
	double InCycFac_6, CInCycFac_6, InCycFac_7, CInCycFac_7, InCycFac_neg_6, CInCycFac_neg_6, InCycFac_neg_7, CInCycFac_neg_7;

	int kon, Ckon;
	int trig, Ctrig;
	int flagdmg, Cflagdmg;
	int flagdmg_Hardening, Cflagdmg_Hardening;
	int flagdmg_Hardening_strength, Cflagdmg_Hardening_strength;
	int flagstop, Cflagstop;
	int flagcap, Cflagcap;

	double dpeakmax, Cdpeakmax;
	double dpeakmax_inner, Cdpeakmax_inner;
	double dpeakmax_inner_inner, Cdpeakmax_inner_inner;
	double dpealmax_bench, Cdpealmax_bench;
	double ffmax_inner_inner, Cffmax_inner_inner;
	double fouterNP_max, CfouterNP_max;
	double ffmin_inner, Cffmin_inner;
	double ffmin_inner_inner, Cffmin_inner_inner;
	double fouterPN_min, CfouterPN_min;
	double dpeakmin, Cdpeakmin;
	double dpeakmin_inner, Cdpeakmin_inner;
	double dpeakmin_inner_inner, Cdpeakmin_inner_inner;
	double dpeakmin_bench, Cdpeakmin_bench;
	double ffmin, Cffmin;
	double ffmax, Cffmax;
	double ffmax_inner, Cffmax_inner;

	double spos, Cspos;
	double sneg, Csneg;
	double foffpos, Cfoffpos;
	double foffneg, Cfoffneg;
	double dmgSpos, CdmgSpos;
	double dmgSneg, CdmgSneg;
	double alpha, Calpha;
	double Et, CEt;
	double Ec, CEc;
	double Eposnorm, CEposnorm;
	double Enegnorm, CEnegnorm;

	double dypos, dyneg;
	double drespos, dresneg, Cdrespos, Cdresneg;
	double dUpos, dUneg;
	double dppos, dpneg, Cdppos, Cdpneg;
	double dpcpos, dpcneg;
	double ahardpos, ahardneg;
	double adegpos, adegneg;
	double Ed0pos, Ed0neg;
	double Ed1pos, Ed1neg;

	int getStateFlag(void);
	void defineBackbone(void);
	void checkEnvelope(void);
	void define_peak(void);
	void update_damage(void);
	void update_damage_hardeingin(void);
	void Total_area(void);
	//-------------------My Variable--------------------------------

								 //Trial Variable
	double TstrainRate;
	double CstrainRate;
	double TstrainMax;		// Maximum strain
	double CstrainMax;
	double TstrainMin;
	double CstrainMin;
	int TstateFlag;
	int CstateFlag;
	double deltaD;
	double CdeltaD;
	double d12;
	double Cd12;
	double f12;
	double Cf12;
	double d_12;
	double Cd_12;
	double f_12;
	double Cf_12;
	double d12neg;
	double Cd12neg;
	double f12neg;
	double Cf12neg;
	double d_12neg;
	double Cd_12neg;
	double f_12neg;
	double Cf_12neg;
	double d2;
	double Cd2;
	double f2;
	double Cf2;
	double d_2;
	double Cd_2;
	double f_2;
	double Cf_2;
	double d2neg;
	double Cd2neg;
	double f2neg;
	double Cf2neg;
	double d_2neg;
	double Cd_2neg;
	double f_2neg;
	double Cf_2neg;
	double d3;
	double Cd3;
	double f3;
	double Cf3;
	double d_3;
	double Cd_3;
	double f_3;
	double Cf_3;
	double d3neg;
	double Cd3neg;
	double f3neg;
	double Cf3neg;
	double d_3neg;
	double Cd_3neg;
	double f_3neg;
	double Cf_3neg;
	double Tdu;
	double CTdu;
	double dmgCounter;
	double CdmgCounter;

	double Benchmark56_up, CBenchmark56_up;
	double Benchmark56_down, CBenchmark56_down;
	double Benchmark_neg_56_up, CBenchmark_neg_56_up;
	double Benchmark_neg_56_down, CBenchmark_neg_56_down;
	double Benchmark67_up, CBenchmark67_up;
	double Benchmark67_down, CBenchmark67_down;
	double Benchmark_neg_67_up, CBenchmark_neg_67_up;
	double Benchmark_neg_67_down, CBenchmark_neg_67_down;

	double spos_pos, Cspos_pos;
	double sneg_pos, Csneg_pos;
	double Et_pos, CEt_pos;
	double Ec_pos, CEc_pos;
	double spos_neg, Cspos_neg;
	double sneg_neg, Csneg_neg;
	double Et_neg, CEt_neg;
	double Ec_neg, CEc_neg;
	double Epos_pos, CEpos_pos;
	double Eneg_neg, CEneg_neg;
	double T_area, CT_area;

	//Hardening damage model
	double Et_hard, CEt_hard;
	double Ec_hard, CEc_hard;
	double spos_hard, Cspos_hard;
	double sneg_hard, Csneg_hard;
	double spos_pos_hard, Cspos_pos_hard;
	double sneg_pos_hard, Csneg_pos_hard;
	double Et_pos_hard, CEt_pos_hard;
	double Ec_pos_hard, CEc_pos_hard;
	double spos_neg_hard, Cspos_neg_hard;
	double sneg_neg_hard, Csneg_neg_hard;
	double Et_neg_hard, CEt_neg_hard;
	double Ec_neg_hard, CEc_neg_hard;
	double Epos_pos_hard, CEpos_pos_hard;
	double Eneg_neg_hard, CEneg_neg_hard;
	double Ed0neg_hard, Ed0pos_hard;
	double Ed1pos_hard, Ed1neg_hard;
	double Epos_hard, Eneg_hard, CEpos_hard, CEneg_hard;
	double Eposnorm_hard, CEposnorm_hard;
	double Enegnorm_hard, CEnegnorm_hard;
	double delta_pos_hard, Cdelta_pos_hard;
	double delta_pos_max_hard, Cdelta_pos_max_hard;
	double delta_neg_hard, Cdelta_neg_hard;
	double delta_neg_max_hard, Cdelta_neg_max_hard;
	double alpha_pos, Calpha_pos;
	double alpha_neg, Calpha_neg;

	//Backbone 
	double slope_pos, Cslope_pos;
	double Intcpt_slope_pos, CIntcpt_slope_pos;
	double BenMark;
	double R_dcapneg;
	double CR_dcapneg;
	double R_dresneg;
	double CR_dresneg;
	double R_drespos;
	double CR_drespos;
	double CIntcpt_deg_pos;
	double slope_neg, Cslope_neg;
	double Intcpt_slope_neg, CIntcpt_slope_neg;
	double Intcpt_deg_pos;
	double Intcpt_deg_neg;
	double CIntcpt_deg_neg;

	double Intcpt_res_pos;
	double CIntcpt_res_pos;
	double Intcpt_res_neg;
	double CIntcpt_res_neg;
	double Intcpt_Xaxis_pos;
	double CIntcpt_Xaxis_pos;
	double Intcpt_Xaxis_neg;
	double CIntcpt_Xaxis_neg;

	double R_dypos, CR_dypos;
	double R_fypos, CR_fypos;
	double R_fcappos, CR_fcappos;
	double R_fyneg, CR_fyneg;
	double R_dyneg, CR_dyneg;
	double R_fcapneg, CR_fcapneg;
	double R_frespos, CR_frespos;
	double R_dppos, CR_dppos;
	double R_fresneg, CR_fresneg;
	double R_Kdegneg, CR_Kdegneg;
	double R_Kdegpos, CR_Kdegpos;
	double R_dcappos, CR_dcappos;
	double pnote;
	int countGlobalEnv;
	// commitstate variable
	double ratio;
	double Cratio;

  //double a, b, c, ee, ff;

  //int commitCalledOnce;		
};

#endif
