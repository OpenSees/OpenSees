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

// Written by: Rasool Ghorbani, University of Texas at San Antonio
// Created: Jun 2023
// Revision: A

// Description: This file contains the class implementation for GMG_CMAC2D Model. 
// GMG_CMAC2D model is intended for use within a two-dimensional lumped-plasticity 
// framework and aimed at capturing the lateral cyclic behavior of non-retrofitted 
// and retrofitted reinforced concrete columns subjected to seismic motions, up to 
// complete loss of lateral strength.

// The material is based on newly developed computational tools for decision-oriented reinforced concrete column simulations 
// [Ghorbani2022] Ghorbani, R. (2022). "Computational Framework for Decision-Oriented Reinforced Concrete Column Simulation Capabilities". PhD Dissertation, The University of Texas at San Antonio. <https://www.proquest.com/docview/2702490424?pq-origsite=gscholar&fromopenview=true> 

#ifndef GMG_CMAC2D_h
#define GMG_CMAC2D_h

#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <array>
#include <vector>
class Matrix;
//class vector;


class GMG_CMAC2D : public UniaxialMaterial
{
public:
	GMG_CMAC2D(int tag,
		double H, double B, double C_C,
		double L, double La,
		double fc,
		double fyL, double dbL, double nL_EW, double nL_NS,
		double fyT, double dbT, double n_leg, double S,
		double lb,
		double P_Axial, 
		int PVM_Flag,
		int number_of_averaged_elements_Elastic,
		int number_of_averaged_elements_Hardening,
		int Retrofit_flag,
		double fje,
		double tj,
		double bj);

	GMG_CMAC2D();
	~GMG_CMAC2D();

	const char *getClassType(void) const { return "GMG_CMAC2D"; };

	int setTrialStrain(double strain, double strainRate = 0.0);
	int setTrialStrain(const Vector strain_from_element);
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
	void splineparam_Shr(double MtoRref, double dpeakpos, double dcappose, double dpeakneg, double dcapneg);

	void spline_curve(double dirtag, double P0x, double P0y, double P4x, double P4y, double Krel, double Kun, double E, double x);
	void spline_curve_Shr(double dirtag, double P0x, double P0y, double P4x, double P4y, double Krel, double Kun, double E, double x);

	void P_M_yield(double Axial_Force);
	void P_M_capping(double Axial_Force);
	void P_M_splice(double Axial_Force);
	void Shear_strength(double Axial_Force);

	void P_M_yield_Opt(double Axial_Force);
	void P_M_capping_Opt(double Axial_Force);
	void P_M_splice_Opt(double Axial_Force);
	void Shear_strength_Opt(double Axial_Force);

	void DamageOutput(void);
	void DamageOutput_FailureType(void);
	void DamageOutput_Hardening(void);
	void DamageOutput_PostCapping(void);
	void DamageOutput_PostResidual(void);
	// Fixed Input Material Variables
	//--- MONOTONIC Backbone Rotation Spring (Flexure Critical)
	double Kepos_Rot, Keneg_Rot, fypos_Rot, fyneg_Rot, fcappos_Rot, fcapneg_Rot, dcappos_Rot, dcapneg_Rot;
	double Kdegpos_Rot, Kdegneg_rot, frespos_Rot, fresneg_Rot, delUpos_Rot, delUneg_Rot;

	//--- MONOTONIC Backbone Rotation Spring (Splice Critical)
	double fypos_splice_Rot, fyneg_splice_Rot, fypos_splice_deg_Rot, fyneg_splice_deg_Rot, fcappos_splice_Rot, fcapneg_splice_Rot, dcappos_splice_Rot, dcapneg_splice_Rot;
	double Kdegpos_splice_Rot, Kdegneg_splice_rot, frespos_splice_Rot, fresneg_splice_Rot, delUpos_splice_Rot, delUneg_splice_Rot;


	// Automatic version
	double B;
	double H;
	double C_C;
	double La;
	double L;


	double fc;

	double dbL;     // Longitudinal reinf. diameter
	double nL_EW;      // # longitudinal bar in EW direction 
	double nL_NS;      // # longitudinal bar in NS direction
	double fyL;

	double dbT;     // Transverse reinf. diameter
	double n_leg;        // # trans. reinf. leg
	double S;       // spacing btween trans. reinf
	double fyT;

	double lb;
	double ld;

	double P_Axial;
	double P_Axial_sum;
	double P_Axial_average;
	int PVM_Flag;
	int number_of_averaged_elements_Elastic;
	int number_of_averaged_elements_Hardening;

	int Retrofit_flag;
	double fje;
	double N_Layer;
	double tj;
	double bj;

	Vector P_Axial_Vector;
	Vector CP_Axial_Vector;


	//--- Section Properties
	double Ckd; // nuetral axis

	//--- MONOTONIC Backbone Rotation Spring (Flexure-Shear Critical)
	double dcappos_FlexShr_Rot, dcapneg_FlexShr_Rot;
	double Kdegpos_FlexShr_Rot, Kdegneg_FlexShr_rot, frespos_FlexShr_Rot, fresneg_FlexShr_Rot, delUpos_FlexShr_Rot, delUneg_FlexShr_Rot;

	//Spline-Curve properties (Flexure Critical)
	double alpha_Er_Rot_Post_Yielding, beta_Er_Rot_Post_Yielding;
	double alpha_Er_Rot_Post_Capping, beta_Er_Rot_Post_Capping;
	double Er_Rot_Post_Yielding, Er_Rot_Post_Capping;
	double Kun_Rot_Post_Yielding, Kun_Rot_Post_Capping;
	double Kr_Rot_Post_Yielding, Kr_Rot_Post_Capping;

	//Spline-Curve properties (Splice Critical)
	double alpha_Er_Rot_Post_Capping_splice;
	double beta_Er_Rot_Post_Capping_splice;
	double Er_Rot_Post_Capping_splice;
	double Kun_Rot_Post_Capping_splice;
	double Kr_Rot_Post_Capping_splice;

	//Spline-Curve properties (Flexure-Shear Critical)
	double alpha_Er_Rot_Post_Capping_FlexShr;
	double beta_Er_Rot_Post_Capping_FlexShr;
	double Er_Rot_Post_Capping_FlexShr;
	double Kun_Rot_Post_Capping_FlexShr;
	double Kr_Rot_Post_Capping_FlexShr;

	//--- DAMAGE Properties (Flexure Critical)
	double C1_Rot, C2_Rot, C3_Rot, solpe_post_yielding_Rot, solpe_post_capping_Rot;
	double delta_ratio_max_hard_Rot, Ref_Energy_Coe_Rot;


	//--- DAMAGE Properties (Splice Critical)
	double C1_splice_Rot, C2_splice_Rot, C3_splice_Rot, solpe_post_capping_splice_Rot;

	//--- DAMAGE Properties (Flexure-Shear Critical)
	double C1_FlexShr_Rot, C2_FlexShr_Rot, C3_FlexShr_Rot, solpe_post_capping_FlexShr_Rot;

	// Fixed Input Material Variables
	//--- MONOTONIC Backbone Shear Spring
	double Kepos_Shr, Keneg_Shr, fypos_Shr, fyneg_Shr, fcappos_Shr, fcapneg_Shr, dcappos_Shr, dcapneg_Shr;
	double Kdegpos_Shr, Kdegneg_Shr, frespos_Shr, fresneg_Shr, delUpos_Shr, delUneg_Shr;
	//Curve properties
	double alpha_Er_Shr_Post_Yielding, beta_Er_Shr_Post_Yielding;
	double alpha_Er_Shr_Post_Capping, beta_Er_Shr_Post_Capping;
	double Er_Shr_Post_Yielding, Er_Shr_Post_Capping;
	//double Kun_Shr_Pos_Post_Yielding, Kun_Shr_Pos_Post_Capping, Kun_Shr_Neg_Post_Yielding, Kun_Shr_Neg_Post_Capping;
	double Kun_Shr_Post_Yielding, Kun_Shr_Post_Capping;
	double Kr_Shr_Pos_Post_Yielding, Kr_Shr_Pos_Post_Capping, Kr_Shr_Neg_Post_Yielding, Kr_Shr_Neg_Post_Capping;
	double Kr_Shr_Post_Yielding, Kr_Shr_Post_Capping;
	//--- DAMAGE Properties
	double C1_Shr, C2_Shr, C3_Shr, solpe_post_yielding_Shr, solpe_post_capping_Shr;
	double delta_ratio_max_hard_Shr, Ref_Energy_Coe_Shr;

	//--- MONOTONIC Backbone Axial Spring
	double Kepos_Axil;
	double Keneg_Axil;



	// State Variables
	double d;     // Deformation
	double f;     // Force
	double F_Rot;
	double F_Shr;
	double ek;    // Tangent Stiffness or slope
	double ek_Rot, Cek_Rot;
	double ek_Shr, Cek_Shr;
	double ek_Axi, Cek_Axi;
	double esel;  // Secant to elastic slope
	double MtoRref, KrR, Krel, KuR, Kun, ER, E;
	double CMtoRref, CKrR, CKrel, CKuR, CKun, CER, CE;
	double MtoRref_Shr, KrR_Shr, Krel_Shr, KuR_Shr, Kun_Shr, ER_Shr, E_Shr;
	double CMtoRref_Shr, CKrR_Shr, CKrel_Shr, CKuR_Shr, CKun_Shr, CER_Shr, CE_Shr;
	double x;
	int dirtag;
	double fspl, fderspl, dderspl;
	double din, fin;
	double fouterPN, fouterNP;
	double fouterPN_Shr, fouterNP_Shr;
	double ekouterPN, ekouterNP;
	double ekouterPN_Shr, ekouterNP_Shr;
	double finner, ekinner;
	double finner_Shr, ekinner_Shr;
	double s, cwct, cwcc, Epos, Eneg;
	double Epos_Shr, Eneg_Shr;
	//double Energy_Coe;

	// Ratio of VyE/Vcol for defferentiating between flexure and flexure-shear
	double VyE_to_Vcol_Ratio, CVyE_to_Vcol_Ratio;
	double Vcap_for_ER, CVcap_for_ER;

	// Cyclic parameters
	// Flexure
	double n;
	double Lamda_KsecS;
	double Lamda_KsecG;
	double Lamda_KrelG;
	double Lamda_KunG;
	double Lamda_KrelS;
	double Lamda_KunS;

	double Lambda_ElasG;
	// Shear
	double Lamda_KsecS_Shr;
	double Lamda_KsecG_Shr;
	double Lamda_KrelG_Shr;
	double Lamda_KunG_Shr;
	double Lamda_KrelS_Shr;
	double Lamda_KunS_Shr;

	// Trial and Committeed State Variables
	double d_Axial, d_Shr, d_Rot;
	double Cstrain;         // Committed Strain
	double Cstrain_Rot;
	double Cstrain_Shr;
	double Cstrain_Axial;
	double Cstress;         // Committed Stress
	double Cstress_Rot;
	double Cstress_Shr;
	double Cstress_Axial, F_Axial;
	double Cek;				// Committed Tangent
	double InCycFac, CInCycFac; // In-cycle Factor
	double InCycFac_Rot, CInCycFac_Rot; // In-cycle Factor
	double InCycFac_Shr, CInCycFac_Shr; // In-cycle Factor

	double InCycFac_6, CInCycFac_6, InCycFac_7, CInCycFac_7, InCycFac_neg_6, CInCycFac_neg_6, InCycFac_neg_7, CInCycFac_neg_7;
	double InCycFac_6_Shr, CInCycFac_6_Shr, InCycFac_7_Shr, CInCycFac_7_Shr, InCycFac_neg_6_Shr, CInCycFac_neg_6_Shr, InCycFac_neg_7_Shr, CInCycFac_neg_7_Shr;

	int kon, Ckon;
	int trig, Ctrig;

	int flagdmg, Cflagdmg;
	int flagdmg_Shr, Cflagdmg_Shr;

	int flagdmg_Hardening_strength, Cflagdmg_Hardening_strength;
	int flagdmg_Hardening_strength_Shr, Cflagdmg_Hardening_strength_Shr;

	int flagdmg_Hardening, Cflagdmg_Hardening;
	int flagdmg_Hardening_Shr, Cflagdmg_Hardening_Shr;

	int flag_entering_hardening_Flex_Rot, Cflag_entering_hardening_Flex_Rot;
	int flag_entering_hardening_Flex_pos_Rot, Cflag_entering_hardening_Flex_pos_Rot;
	int flag_entering_hardening_Flex_neg_Rot, Cflag_entering_hardening_Flex_neg_Rot;
	int flag_entering_hardening_Shr_pos, Cflag_entering_hardening_Shr_pos;
	int flag_entering_hardening_Shr_neg, Cflag_entering_hardening_Shr_neg;
	int flag_entering_hardening_Pos_Flex_Rot, Cflag_entering_hardening_Pos_Flex_Rot;
	int flag_entering_hardening_Neg_Flex_Rot, Cflag_entering_hardening_Neg_Flex_Rot;
	int flag_entering_hardening_Splice_Rot, Cflag_entering_hardening_Splice_Rot;
	int flag_entering_hardening_Splice_pos_Rot, Cflag_entering_hardening_Splice_pos_Rot;
	int flag_entering_hardening_Splice_neg_Rot, Cflag_entering_hardening_Splice_neg_Rot;
	int flag_entering_hardening_Shr, Cflag_entering_hardening_Shr;
	int flag_entering_residual_Flex_Rot, Cflag_entering_residual_Flex_Rot;
	int flag_entering_residual_Shr, Cflag_entering_residual_Shr;

	int flag_FlexShr_Filure_Rot, Cflag_FlexShr_Filure_Rot;
	int flag_Splice_Filure_Rot, Cflag_Splice_Filure_Rot;
	int flag_FlexSplice_Filure_Rot, Cflag_FlexSplice_Filure_Rot;
	int flag_Flx_Filure_Rot, Cflag_Flx_Filure_Rot;
	int flag_Filure_Shr, Cflag_Filure_Shr;

	int flagstop, Cflagstop;
	int flagcap, Cflagcap;

	double dpeakmax, Cdpeakmax;
	double dpeakmax_Shr, Cdpeakmax_Shr;

	double dpeakmax_bench, Cdpeakmax_bench;
	double dpeakmax_bench_Shr, Cdpeakmax_bench_Shr;


	double d_Bench_Rot, Cd_Bench_Rot;
	double d_Bench_Shr, Cd_Bench_Shr;

	double f_Bench_Rot, Cf_Bench_Rot;
	double f_Bench_Shr, Cf_Bench_Shr;

	int flagFlur_Shr, CflagFlur_Shr;
	int flagFlur_Rot, CflagFlur_Rot;

	double dpeakmax_inner, Cdpeakmax_inner;
	double dpeakmax_inner_Shr, Cdpeakmax_inner_Shr;

	double dpeakmax_inner_inner, Cdpeakmax_inner_inner;
	double dpeakmax_inner_inner_Shr, Cdpeakmax_inner_inner_Shr;

	double ffmax_inner_inner, Cffmax_inner_inner;
	double ffmax_inner_inner_Shr, Cffmax_inner_inner_Shr;

	double fouterNP_max, CfouterNP_max;
	double fouterNP_max_Shr, CfouterNP_max_Shr;

	double ffmin_inner, Cffmin_inner;
	double ffmin_inner_Shr, Cffmin_inner_Shr;

	double ffmin_inner_inner, Cffmin_inner_inner;
	double ffmin_inner_inner_Shr, Cffmin_inner_inner_Shr;

	double fouterPN_min, CfouterPN_min;
	double fouterPN_min_Shr, CfouterPN_min_Shr;

	double dpeakmin, Cdpeakmin;
	double dpeakmin_Shr, Cdpeakmin_Shr;

	double dpeakmin_bench, Cdpeakmin_bench;
	double dpeakmin_bench_Shr, Cdpeakmin_bench_Shr;

	double dpeakmin_inner, Cdpeakmin_inner;
	double dpeakmin_inner_Shr, Cdpeakmin_inner_Shr;

	double dpeakmin_inner_inner, Cdpeakmin_inner_inner;
	double dpeakmin_inner_inner_Shr, Cdpeakmin_inner_inner_Shr;

	double ffmin, Cffmin;
	double ffmin_Shr, Cffmin_Shr;

	double ffmax, Cffmax;
	double ffmax_Shr, Cffmax_Shr;

	double ffmax_inner, Cffmax_inner;
	double ffmax_inner_Shr, Cffmax_inner_Shr;

	double spos, Cspos;
	double spos_Shr, Cspos_Shr;

	double sneg, Csneg;
	double sneg_Shr, Csneg_Shr;

	double foffpos, Cfoffpos;
	double foffneg, Cfoffneg;
	double dmgSpos, CdmgSpos;
	double dmgSpos_Shr, CdmgSpos_Shr;

	double dmgSneg, CdmgSneg;
	double dmgSneg_Shr, CdmgSneg_Shr;

	double Et, CEt;
	double Et_Shr, CEt_Shr;

	double Ec, CEc;
	double Ec_Shr, CEc_Shr;

	double Eposnorm, CEposnorm;
	double Eposnorm_Shr, CEposnorm_Shr;

	double Enegnorm, CEnegnorm;
	double Enegnorm_Shr, CEnegnorm_Shr;


	// Static Variables
	double Vcol, CVcol;
	double dypos, dyneg;
	double dypos_splice, dyneg_splice;
	double dypos_Shr, dyneg_Shr;
	double drespos, dresneg;
	double drespos_Shr, dresneg_Shr;
	double dUpos, dUneg;
	double dUpos_Shr, dUneg_Shr;

	double fs, lb_deg, fs_deg;

	double dppos, dpneg, Cdppos, Cdpneg;
	double dppos_Shr, dpneg_Shr, Cdppos_Shr, Cdpneg_Shr;

	double dpcpos, dpcneg;
	double dpcpos_Shr, dpcneg_Shr;

	double ahardpos, ahardneg;
	double ahardpos_Shr, ahardneg_Shr;

	double adegpos, adegneg;
	double adegpos_Shr, adegneg_Shr;

	double Ed0pos, Ed0neg;
	double Ed0pos_Shr, Ed0neg_Shr;

	double Ed1pos, Ed1neg;
	double Ed1pos_Shr, Ed1neg_Shr;


	//--------------------Adding my function------------------------
	int getStateFlag(void);
	int getStateFlag_Shr(void);

	void defineBackbone(void);
	void defineBackbone_Shr(void);

	void Calibration(void);
	void Calibration_Shr(void);
	void Needed_Inpt_to_get_failure_mode(void);

	void checkEnvelope(void);
	void checkEnvelope_Shr(void);

	void define_peak(void);
	void define_peak_Shr(void);

	void update_damage(void);
	void update_damage_Shr(void);

	void update_damage_hardeingin(void);
	void update_damage_hardeingin_Shr(void);

	//void Total_area(void);
	double GMG_CMAC2D::getTangent_Multi(void);
	double GMG_CMAC2D::getStress_Multi(void);
	double GMG_CMAC2D::getInitialTangent_Multi(void);
	double GMG_CMAC2D::getStrain_Multi(void);
	double GMG_CMAC2D::getDampTangent_Multi(void);
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
	int TstateFlag_Shr;
	int CstateFlag_Shr;
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
	double Tdu_Rot, CTdu_Rot;
	double Tdu_Shr, CTdu_Shr;
	double dmgCounter;
	double CdmgCounter;

	double dmgCounter_Shr;
	double CdmgCounter_Shr;

	double Benchmark56, CBenchmark56;
	double Benchmark56_Shr, CBenchmark56_Shr;

	double Benchmark56_up, CBenchmark56_up;
	double Benchmark56_up_Shr, CBenchmark56_up_Shr;

	double Benchmark56_down, CBenchmark56_down;
	double Benchmark56_down_Shr, CBenchmark56_down_Shr;

	double Benchmark_neg_56_up, CBenchmark_neg_56_up;
	double Benchmark_neg_56_up_Shr, CBenchmark_neg_56_up_Shr;

	double Benchmark_neg_56_down, CBenchmark_neg_56_down;
	double Benchmark_neg_56_down_Shr, CBenchmark_neg_56_down_Shr;

	double Benchmark67_up, CBenchmark67_up;
	double Benchmark67_up_Shr, CBenchmark67_up_Shr;

	double Benchmark67_down, CBenchmark67_down;
	double Benchmark67_down_Shr, CBenchmark67_down_Shr;

	double Benchmark_neg_67_up, CBenchmark_neg_67_up;
	double Benchmark_neg_67_up_Shr, CBenchmark_neg_67_up_Shr;

	double Benchmark_neg_67_down, CBenchmark_neg_67_down;
	double Benchmark_neg_67_down_Shr, CBenchmark_neg_67_down_Shr;

	double spos_pos, Cspos_pos;
	double spos_pos_Shr, Cspos_pos_Shr;

	double sneg_pos, Csneg_pos;
	double sneg_pos_Shr, Csneg_pos_Shr;

	double Et_pos, CEt_pos;
	double Et_pos_Shr, CEt_pos_Shr;

	double Ec_pos, CEc_pos;
	double Ec_pos_Shr, CEc_pos_Shr;

	double spos_neg, Cspos_neg;
	double spos_neg_Shr, Cspos_neg_Shr;

	double sneg_neg, Csneg_neg;
	double sneg_neg_Shr, Csneg_neg_Shr;

	double Et_neg, CEt_neg;
	double Et_neg_Shr, CEt_neg_Shr;

	double Ec_neg, CEc_neg;
	double Ec_neg_Shr, CEc_neg_Shr;

	double Epos_pos, CEpos_pos;
	double Epos_pos_Shr, CEpos_pos_Shr;

	double Eneg_neg, CEneg_neg;
	double Eneg_neg_Shr, CEneg_neg_Shr;

	double T_area, CT_area;
	double T_area_Shr, CT_area_Shr;

	//Hardening damage model Flexure
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
	double  Epos_hard, Eneg_hard, CEpos_hard, CEneg_hard;
	double Eposnorm_hard, CEposnorm_hard;
	double Enegnorm_hard, CEnegnorm_hard;
	double delta_pos_hard, Cdelta_pos_hard;
	double delta_pos_max_hard, Cdelta_pos_max_hard;
	double delta_neg_hard, Cdelta_neg_hard;
	double delta_neg_max_hard, Cdelta_neg_max_hard;
	double alpha_pos, Calpha_pos;
	double alpha_neg, Calpha_neg;

	//Hardening damage model Shear
	double Et_hard_Shr, CEt_hard_Shr;
	double Ec_hard_Shr, CEc_hard_Shr;
	double spos_hard_Shr, Cspos_hard_Shr;
	double sneg_hard_Shr, Csneg_hard_Shr;
	double spos_pos_hard_Shr, Cspos_pos_hard_Shr;
	double sneg_pos_hard_Shr, Csneg_pos_hard_Shr;
	double Et_pos_hard_Shr, CEt_pos_hard_Shr;
	double Ec_pos_hard_Shr, CEc_pos_hard_Shr;
	double spos_neg_hard_Shr, Cspos_neg_hard_Shr;
	double sneg_neg_hard_Shr, Csneg_neg_hard_Shr;
	double Et_neg_hard_Shr, CEt_neg_hard_Shr;
	double Ec_neg_hard_Shr, CEc_neg_hard_Shr;
	double Epos_pos_hard_Shr, CEpos_pos_hard_Shr;
	double Eneg_neg_hard_Shr, CEneg_neg_hard_Shr;
	double Ed0neg_hard_Shr, Ed0pos_hard_Shr;
	double Ed1pos_hard_Shr, Ed1neg_hard_Shr;
	double Epos_hard_Shr, Eneg_hard_Shr, CEpos_hard_Shr, CEneg_hard_Shr;
	double Eposnorm_hard_Shr, CEposnorm_hard_Shr;
	double Enegnorm_hard_Shr, CEnegnorm_hard_Shr;
	double delta_pos_hard_Shr, Cdelta_pos_hard_Shr;
	double delta_pos_max_hard_Shr, Cdelta_pos_max_hard_Shr;
	double delta_neg_hard_Shr, Cdelta_neg_hard_Shr;
	double delta_neg_max_hard_Shr, Cdelta_neg_max_hard_Shr;
	double alpha_pos_Shr, Calpha_pos_Shr;
	double alpha_neg_Shr, Calpha_neg_Shr;

	//Backbone 
	double slope_pos;
	double slope_pos_Shr;

	double Cslope_pos;
	double Cslope_pos_Shr;

	double Intcpt_slope_pos;
	double Intcpt_slope_pos_Shr;

	double CIntcpt_slope_pos;
	double CIntcpt_slope_pos_Shr;

	double BenMark, BenMark_Shr;
	double R_dcapneg;
	double R_dcapneg_Shr;

	double CR_dcapneg;
	double CR_dcapneg_Shr;

	double CR_dcappos;
	double CR_dcappos_Shr;

	double R_dresneg;
	double R_dresneg_Shr;

	double CR_dresneg;
	double CR_dresneg_Shr;

	double R_drespos;
	double R_drespos_Shr;

	double CR_drespos;
	double CR_drespos_Shr;

	double CIntcpt_deg_pos;
	double CIntcpt_deg_pos_Shr;

	double slope_neg;
	double slope_neg_Shr;

	double Cslope_neg;
	double Cslope_neg_Shr;

	double Intcpt_slope_neg;
	double Intcpt_slope_neg_Shr;

	double CIntcpt_slope_neg;
	double CIntcpt_slope_neg_Shr;

	double Intcpt_deg_pos;
	double Intcpt_deg_pos_Shr;

	double Intcpt_deg_neg;
	double Intcpt_deg_neg_Shr;

	double Intcpt_res_pos, CIntcpt_res_pos;
	double Intcpt_res_pos_Shr, CIntcpt_res_pos_Shr;

	double Intcpt_res_neg, CIntcpt_res_neg;
	double Intcpt_res_neg_Shr, CIntcpt_res_neg_Shr;

	double Intcpt_Xaxis_pos, CIntcpt_Xaxis_pos;
	double Intcpt_Xaxis_neg, CIntcpt_Xaxis_neg;

	double Intcpt_Xaxis_pos_Shr, CIntcpt_Xaxis_pos_Shr;
	double Intcpt_Xaxis_neg_Shr, CIntcpt_Xaxis_neg_Shr;

	double CIntcpt_deg_neg;
	double CIntcpt_deg_neg_Shr;


	double R_dypos;
	double R_dypos_Shr;

	double CR_dypos;
	double CR_dypos_Shr;

	double R_fypos;
	double R_fypos_Shr;

	double CR_fypos;
	double CR_fypos_Shr;

	double R_fcappos;
	double R_fcappos_Shr;

	double CR_fcappos;
	double CR_fcappos_Shr;

	double R_fyneg;
	double R_fyneg_Shr;

	double CR_fyneg;
	double CR_fyneg_Shr;

	double R_dyneg;
	double R_dyneg_Shr;

	double CR_dyneg;
	double CR_dyneg_Shr;

	double R_fcapneg;
	double R_fcapneg_Shr;

	double CR_fcapneg;
	double CR_fcapneg_Shr;

	double R_frespos;
	double R_frespos_Shr;

	double CR_frespos;
	double CR_frespos_Shr;

	double R_delUpos;
	double R_delUneg;

	double R_delUpos_Shr;
	double R_delUneg_Shr;

	double CR_delUpos;
	double CR_delUneg;

	double CR_delUpos_Shr;
	double CR_delUneg_Shr;

	double R_dppos;
	double R_fresneg;
	double R_fresneg_Shr;

	double CR_fresneg;
	double CR_fresneg_Shr;

	double R_Kdegpos;
	double R_Kdegpos_Shr;

	double CR_Kdegpos;
	double CR_Kdegpos_Shr;

	double R_Kdegneg;
	double R_Kdegneg_Shr;

	double CR_Kdegneg;
	double CR_Kdegneg_Shr;

	double R_dcappos;
	double R_dcappos_Shr;

	double pnote;
	int countGlobalEnv;
	// commitstate variable
	double ratio;
	double Cratio;

	double ratio_Shr;
	double Cratio_Shr;

	double a, b, c, ee, ff;



	// ++++++++++++++++++ Prepairing the material to act as a Multi axial material ++++++++++++
	Vector Stress_vector;
	Vector Strain_vector;
	Vector Tang_vector;
	Vector Damp_vector;

	Matrix Force_section_yield;
	Matrix Moment_Section_yield;
	Matrix Force_section_capping;
	Matrix Moment_Section_capping;
	Matrix Force_section_splice;
	Matrix Moment_Section_splice;

	static Matrix Damage_Data;
	static int Counter_MACE;
};

#endif