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
#include <math.h>
#include <string.h>
//#include <boost/math/distribution/beta.hpp>

#include <elementAPI.h>
#include <GMG_CyclicReinforcedConcrete.h>
#include <Vector.h>
#include <Channel.h>
#include <MaterialResponse.h>

#include <OPS_Globals.h>

// Static Variables
static int numGMG_CyclicReinforcedConcreteMaterials = 0;

#ifndef fmax
#define max(a,b) (((a)>(b)) ? (a) : (b))
#endif

#ifndef fmin
#define min(a,b) (((a)<(b)) ? (a) : (b))
#endif

void *
OPS_GMG_CyclicReinforcedConcrete()
{
	if (numGMG_CyclicReinforcedConcreteMaterials == 0) {
		numGMG_CyclicReinforcedConcreteMaterials++;
		opserr << "GMG_CyclicReinforcedConcrete Material Model\n";
		opserr << "Written by R. Ghorbani UTSA Copyright 2022\n";
	}

	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial *theMaterial = 0;
	int argc = OPS_GetNumRemainingInputArgs();

	if (!(argc == 32)) {
		opserr << "WARNING GMG_CyclicReinforcedConcreteMaterial -- insufficient arguments\n";
		opserr << "For direct input, the material needs:\n";
		opserr << "UniaxialMaterial GMG_CyclicReinforcedConcrete matTag? Kepos? Keneg?\n";
		opserr << "fypos? fyneg? fcappos? fcapneg? dcappos? dcapneg?\n";
		opserr << "Kdegpos? Kdegneg? frespos? fresneg? delUpos? delUneg?\n";
		opserr << "alpha_Er_Hardening? beta_Er_Hardening?\n";
		opserr << "alpha_Er_Post_Capping? beta_Er_Post_Capping?\n";
		opserr << "ErMax_Hardening? ErMax_Post_Capping?\n";
		opserr << "alpha_Kun_Hardening? alpha_Kun_Post_Capping? beta_Krel_Hardening? beta_Krel_Post_Capping?\n";
		opserr << "delta_ratio_max_hard? Ref_Energy_Coe? C1? C2? C3? solpe_damage_Hardening? solpe_damage_post_cappin?\n";
		return 0;
	}


	int iTagData[1];
	double dMonoData[24];
	double dDmgData[7];
	int numData;

	numData = 1;
	if (OPS_GetIntInput(&numData, iTagData) != 0) {
		opserr << "WARNING GMG_CyclicReinforcedConcreteMaterial -- invalid uniaxialMaterial matTag\n";
		return 0;
	}
	numData = 24;
	if (OPS_GetDoubleInput(&numData, dMonoData)) {
		opserr << "WARNING GMG_CyclicReinforcedConcreteMaterial -- invalid uniaxialMaterial Backbone Properties\n";
		opserr << "For direct input, the material needs:\n";
		opserr << "UniaxialMaterial GMG_CyclicReinforcedConcrete matTag? Kepos? Keneg?\n";
		opserr << "fypos? fyneg? fcappos? fcapneg? dcappos? dcapneg?\n";
		opserr << "Kdegpos? Kdegneg? frespos? fresneg? delUpos? delUneg?\n";
		opserr << "alpha_Er_Hardening? beta_Er_Hardening?\n";
		opserr << "alpha_Er_Post_Capping? beta_Er_Post_Capping?\n";
		opserr << "ErMax_Hardening? ErMax_Post_Capping?\n";
		opserr << "alpha_Kun_Hardening? alpha_Kun_Post_Capping? beta_Krel_Hardening? beta_Krel_Post_Capping?\n";
		return 0;
	}

	numData = 7;
	if (OPS_GetDoubleInput(&numData, dDmgData)) {
		opserr << "WARNING GMG_CyclicReinforcedConcreteMaterial -- invalid uniaxialMaterial Damage Properties\n";
		opserr << "For direct input, the material needs:\n";
		opserr << "delta_ratio_max_hard? Ref_Energy_Coe? C1? C2? C3? solpe_damage_Hardening? solpe_damage_post_cappin?\n";
		return 0;
	}
	// Parsing was successful, allocate the material
	theMaterial = new GMG_CyclicReinforcedConcrete(iTagData[0], dMonoData[0], dMonoData[1], dMonoData[2],
		dMonoData[3], dMonoData[4], dMonoData[5], dMonoData[6], dMonoData[7], dMonoData[8],
		dMonoData[9], dMonoData[10], dMonoData[11], dMonoData[12], dMonoData[13],
		dMonoData[14], dMonoData[15],
		dMonoData[16], dMonoData[17], dMonoData[18], dMonoData[19],
		dMonoData[20], dMonoData[21], dMonoData[22], dMonoData[23],
		dDmgData[0], dDmgData[1], dDmgData[2], dDmgData[3], dDmgData[4], dDmgData[5], dDmgData[6]);
	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type GMG_CyclicReinforcedConcrete Material\n";
		return 0;
	}
	return theMaterial;
}

// Constructors
GMG_CyclicReinforcedConcrete::GMG_CyclicReinforcedConcrete(int tag, double Kepos, double Keneg, double fypos, double fyneg,
	double fcappos, double fcapneg, double dcappos, double dcapneg,
	double Kdegpos, double Kdegneg,
	double frespos, double fresneg,
	double delUpos, double delUneg,
	double alpha_Er_Post_Yielding, double beta_Er_Post_Yielding,
	double alpha_Er_Post_Capping, double beta_Er_Post_Capping,
	double Er_Post_Yielding, double Er_Post_Capping,
	double Kun_Post_Yielding, double Kun_Post_Capping, double Kr_Post_Yielding, double Kr_Post_Capping,
	double delta_ratio_max_hard, double Ref_Energy_Coe,
	double C1, double C2, double C3, double solpe_post_yielding, double solpe_post_capping)
	:UniaxialMaterial(tag, MAT_TAG_GMG_CyclicReinforcedConcrete), Kepos(Kepos), Keneg(Keneg), fypos(fypos), fyneg(fyneg), fcappos(fcappos),
	fcapneg(fcapneg), dcappos(dcappos), dcapneg(dcapneg), Kdegpos(Kdegpos), Kdegneg(Kdegneg), frespos(frespos), fresneg(fresneg),
	delUpos(delUpos), delUneg(delUneg),
	alpha_Er_Post_Yielding(alpha_Er_Post_Yielding), beta_Er_Post_Yielding(beta_Er_Post_Yielding),
	alpha_Er_Post_Capping(alpha_Er_Post_Capping), beta_Er_Post_Capping(beta_Er_Post_Capping),
	Er_Post_Yielding(Er_Post_Yielding), Er_Post_Capping(Er_Post_Capping),
	Kun_Post_Yielding(Kun_Post_Yielding), Kun_Post_Capping(Kun_Post_Capping), Kr_Post_Yielding(Kr_Post_Yielding), Kr_Post_Capping(Kr_Post_Capping),
	delta_ratio_max_hard(delta_ratio_max_hard), Ref_Energy_Coe(Ref_Energy_Coe),
	C1(C1), C2(C2), C3(C3), solpe_post_yielding(solpe_post_yielding), solpe_post_capping(solpe_post_capping)
{
	// Initialize Variables by Calling revertToStart function
	this->revertToStart();
}

GMG_CyclicReinforcedConcrete::GMG_CyclicReinforcedConcrete()
	:UniaxialMaterial(0, MAT_TAG_GMG_CyclicReinforcedConcrete), Kepos(0.0), Keneg(0.0), fypos(0.0), fyneg(0.0), fcappos(0.0),
	fcapneg(0.0), dcappos(0.0), dcapneg(0.0), Kdegpos(0.0), Kdegneg(0.0), frespos(0.0), fresneg(0.0),
	delUpos(0.0), delUneg(0.0),
	alpha_Er_Post_Yielding(0.0), beta_Er_Post_Yielding(0.0),
	alpha_Er_Post_Capping(0.0), beta_Er_Post_Capping(0.0),
	Er_Post_Yielding(0.0), Er_Post_Capping(0.0),
	Kun_Post_Yielding(0.0), Kun_Post_Capping(0.0), Kr_Post_Yielding(0.0), Kr_Post_Capping(0.0),
	delta_ratio_max_hard(0.0), Ref_Energy_Coe(0.0),	 
	C1(0.0), C2(0.0), C3(0.0), solpe_post_yielding(0.0), solpe_post_capping(0.0)

{
	// Initialize Variables by Calling revertToStart function
	this->revertToStart();
	//revertToStart();
}

GMG_CyclicReinforcedConcrete::~GMG_CyclicReinforcedConcrete()
{
	// does nothing
}

int
GMG_CyclicReinforcedConcrete::setTrialStrain(double strain, double strainRate)
{
	// all variables to the last commit state
	this->revertToLastCommit();
	//revertToLastCommit();
	d = strain;
	TstrainRate = strainRate;
	Tdu = d - Cstrain;

	if (Tdu == 0.0 || fabs(Tdu) > 1.0)
		return 0;

	if (d >= TstrainMax)
		TstrainMax = d;
	else if (d < TstrainMin)
		TstrainMin = d;

	if (TstateFlag == 0)
	{
		flagdmg = 0;
		if (d >= 0.0) {
			f = Kepos * d;
			ek = Kepos;
			if (d >= dpeakmax) {
				TstateFlag = 12;
				defineBackbone();
				f = slope_pos * d + Intcpt_slope_pos;
				ek = slope_pos;

			}

		}
		else if (d < 0.0) {
			f = Keneg * d;
			ek = Keneg;
			if (d <= dpeakmin) {
				TstateFlag = -12;
				defineBackbone();
				f = -(slope_neg * fabs(d) - Intcpt_slope_neg);
				ek = slope_neg;

			}
		}
	}

	else
	{
		TstateFlag = getStateFlag();
		switch (TstateFlag)
		{
		case 1:	f = slope_pos * d + Intcpt_slope_pos;
			ek = slope_pos;
			break;
		case -1: f = -(slope_neg * fabs(d) - Intcpt_slope_neg);
			ek = slope_neg;
			break;
		case 12:
			flagdmg = 0;
			flagdmg_Hardening = 1;
			flagdmg_Hardening_strength = 1;

			if (BenMark > 0) {
				ek = slope_pos;
				f = slope_pos * d + Intcpt_slope_pos;
			}
			else if (BenMark < 0) {
				ek = slope_pos;
				f = slope_pos * fabs(d) + Intcpt_slope_pos;
			}
			dpeakmax = d;
			if (flagdmg_Hardening == 1)
				update_damage_hardeingin();

			if (flagdmg_Hardening_strength == 1)
				update_damage();
			break;

		case -12:
			flagdmg = 0;
			flagdmg_Hardening = 1;
			flagdmg_Hardening_strength = 1;
			if (BenMark > 0) {
				ek = slope_neg;
				f = -(slope_neg * fabs(d) - Intcpt_slope_neg);
				d12neg = d;
				f12neg = f;
			}
			else if (BenMark < 0) {
				ek = slope_neg;
				f = -(slope_neg * fabs(d) - Intcpt_slope_neg);
				d_12neg = d;
				f_12neg = f;
			}
			dpeakmin = d;
			if (flagdmg_Hardening == 1)
				update_damage_hardeingin();
			if (flagdmg_Hardening_strength == 1)
				update_damage();
			break;

		case 2:
			//---------Damage will start after entering case 2 or case -2 for the first time---------
			flagdmg = 1;
			flagdmg_Hardening = 0;
			flagdmg_Hardening_strength = 0;
			if (BenMark > 0) {
				ek = Kdegpos;
				f = Kdegpos * fabs(d) + Intcpt_deg_pos;
			}
			else if (BenMark < 0) {
				ek = Kdegpos;
				f = Kdegpos * fabs(d) + Intcpt_deg_pos;
				d_2 = d;
				f_2 = f;
			}
			if (flagdmg == 1)
				update_damage();
			break;

		case -2:
			//---------Damage will start after entering case 2 or case -2 for the first time---------
			flagdmg = 1;
			flagdmg_Hardening = 0;
			flagdmg_Hardening_strength = 0;
			if (BenMark > 0) {
				ek = Kdegneg;
				f = -(Kdegneg*fabs(d) + Intcpt_deg_neg);
				d2neg = d;
				f2neg = f;
			}
			else if (BenMark < 0) {
				ek = Kdegneg;
				f = -(Kdegneg * fabs(d) + Intcpt_deg_neg);
				d_2neg = d;
				f_2neg = f;

			}
			if (flagdmg == 1)
				update_damage();
			break;

		case 3:
			flagdmg = 0;
			flagdmg_Hardening = 0;
			flagdmg_Hardening_strength = 0;
			if (BenMark > 0) {
				ek = -0.0001;
				f = R_frespos;
				d3 = d;
				f3 = f;
			}
			else if (BenMark < 0) {
				ek = -0.0001;
				f = R_frespos;
				d_3 = d;
				f_3 = f;
			}
			break;

		case -3:
			flagdmg = 0;
			flagdmg_Hardening = 0;
			flagdmg_Hardening_strength = 0;
			if (BenMark > 0) {
				ek = -0.0001;
				f = R_fresneg;
				d3neg = d;
				f3neg = f;
			}
			else if (BenMark < 0) {
				ek = -0.0001;
				f = R_fresneg;
				d_3neg = d;
				f_3neg = f;
			}
			break;

		case 30:
			flagdmg = 0;
			flagdmg_Hardening = 0;
			flagdmg_Hardening_strength = 0;
			if (BenMark > 0) {
				ek = Kdegpos;
				f = Kdegpos * fabs(d) + Intcpt_res_pos;
			}
			else if (BenMark < 0) {
				ek = Kdegpos;
				f = Kdegpos * fabs(d) + Intcpt_res_pos;
			}
			break;

		case -30:
			flagdmg = 0;
			flagdmg_Hardening = 0;
			flagdmg_Hardening_strength = 0;
			if (BenMark > 0) {
				ek = Kdegneg;
				f = -(Kdegneg * fabs(d) + Intcpt_res_neg);
			}
			else if (BenMark < 0) {
				ek = Kdegneg;
				f = -(Kdegneg * fabs(d) + Intcpt_res_neg);
			}
			break;


		case 40:
			flagdmg = 0;
			flagdmg_Hardening = 0;
			flagdmg_Hardening_strength = 0;
			if (BenMark > 0) {
				ek = -0.0001;
				f = 0.0;
			}
			else if (BenMark < 0) {
				ek = -0.0001;
				f = 0.0;
			}
			break;

		case -40:

			flagdmg = 0;
			flagdmg_Hardening = 0;
			flagdmg_Hardening_strength = 0;
			if (BenMark > 0) {
				ek = -0.0001;
				f = 0.0;
			}
			else if (BenMark < 0) {
				ek = -0.0001;
				f = 0.0;
			}
			break;

		case 31:
			if (CstateFlag != 31)
				define_peak();
			ek = (ffmax - ffmin) / (dpeakmax - dpeakmin);
			f = ek * (d - dpeakmax) + ffmax;
			checkEnvelope();
			break;

		case 41:
			ek = 0.0001;
			f = 0.0;
			break;


		case 4:

			define_peak();
			dirtag = 1; // Positive to Negative Direction
			MtoRref = (ffmax - ffmin) / (dpeakmax - dpeakmin);
			splineparam(MtoRref, dpeakmax, R_dcappos, dpeakmin, R_dcapneg);
			Krel = pow(MtoRref / ((Kepos + Keneg) / 2), KrR) * (Kepos + Keneg) / 2;
			Kun = fmin(KuR * MtoRref, (Kepos + Keneg) / 2);
			E = ER * (ffmax - ffmin) * ((dpeakmax - dpeakmin) - (ffmax - ffmin) / Kepos);  
			spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d);      
			f = fspl;

			if (flagdmg == 1)
				update_damage();
			if (flagdmg_Hardening == 1)
				update_damage_hardeingin();
			if (flagdmg_Hardening_strength == 1)
				update_damage();
			dpeakmin_inner = d;
			ffmin_inner = f;
			break;

		case 5:
			spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d);
			f = fspl;
			if (flagdmg == 1)
				update_damage();
			if (flagdmg_Hardening == 1)
				update_damage_hardeingin();
			if (flagdmg_Hardening_strength == 1)
				update_damage();
			checkEnvelope();
			break;

		case 6:

			define_peak();
			// Because of the damage model our peaks may be updated so we need to re-define spline parametere again 
			if (CstateFlag != 6) {
				MtoRref = (ffmax - ffmin) / (dpeakmax - dpeakmin);
				splineparam(MtoRref, dpeakmax, R_dcappos, dpeakmin, R_dcapneg);
				Krel = pow(MtoRref / ((Kepos + Keneg) / 2), KrR) * (Kepos + Keneg) / 2;
				Kun = fmin(KuR * MtoRref, (Kepos + Keneg) / 2);
				E = ER * (ffmax - ffmin) * ((dpeakmax - dpeakmin) - (ffmax - ffmin) / Keneg);
			}
			// -----I needed to check if my current point is above or under the new lower limit------
			if (CstateFlag != 6) {
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain); // NP
				Benchmark56_up = fspl;

				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain); // PN
				Benchmark56_down = fspl;
			}
			//---------------------------------------------------------------------------------------
			if ((ffmin_inner <= Benchmark56_up && ffmin_inner >= Benchmark56_down)) {

				if (CstateFlag != 6) {

					x = Cstrain;
					spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
					fouterPN_min = fspl;
					ekouterPN = ek;

					x = Cstrain;
					spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
					fouterNP = fspl;
					ekouterNP = ek;

					InCycFac = (fouterNP - ffmin_inner) / (fouterNP - fouterPN_min);

				}

				x = d;
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
				fouterNP = fspl;
				ekouterNP = ek;
				x = dpeakmin + (d - dpeakmin_inner);
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
				finner = fouterNP - InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)));
				ekinner = (ek - ekouterNP) / ((fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) * (InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) + ekouterNP;
				ek = ekinner;
				f = fmin(finner, fouterNP);
				if (f == fouterNP) {
					ek = ekouterNP;
					TstateFlag = -5;
				}
			}


			if (ffmin_inner < Benchmark56_down) {
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d);
				fouterPN = fspl;
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d);
				fouterNP = fspl;
				ekouterNP = ek;
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmax); // PN
				f = ffmin_inner + ek * (d - dpeakmin_inner);

				f = fmin(f, fouterNP);
				if (f == fouterNP) {
					ek = ekouterNP;
					TstateFlag = -5;
				}

				// we should check if it has reached the lower limit or not, because after it passes the lower limit we should start the coppy-pasting process
				if (f >= fouterPN && Cstress < fouterPN && f < fouterNP) {
					//f = fouterNP;
					dpeakmin_inner = d;
					ffmin_inner = f;
					x = d;
					spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
					fouterPN_min = fspl;
					ekouterPN = ek;

					x = d;
					spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
					fouterNP = fspl;
					ekouterNP = ek;

					InCycFac = (fouterNP - ffmin_inner) / (fouterNP - fouterPN_min);
				}

				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmin); // NP

				if (f >= fouterPN && Cstress > fouterPN && f < fouterNP) {
					x = d;
					spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
					fouterNP = fspl;
					ekouterNP = ek;
					x = dpeakmin + (d - dpeakmin_inner);
					spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
					finner = fouterNP - InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)));
					ekinner = (ek - ekouterNP) / ((fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) * (InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) + ekouterNP;
					ek = ekinner;
					f = fmin(finner, fouterNP);
					if (f == fouterNP) {
						ek = ekouterNP;
						TstateFlag = -5;
					}
				}
			}

			if (ffmin_inner > Benchmark56_up) {
				if (CstateFlag != 6) {
					InCycFac = (ffmin_inner - Benchmark56_up);
				}
				x = d;
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP 
				f = fspl + (InCycFac + (-InCycFac / (dpeakmax - dpeakmin_inner))*(d - dpeakmin_inner));
				ek = ek - (InCycFac / (dpeakmax - dpeakmin_inner));
			}

			if (flagdmg == 1)
				update_damage();
			if (flagdmg_Hardening == 1)
				update_damage_hardeingin();
			if (flagdmg_Hardening_strength == 1)
				update_damage();

			checkEnvelope();
			break;

		case 7:
			define_peak();

			// Because of the damage model our peaks may be updated so we need to re-define spline parametere again 
			if (CstateFlag != 7) {
				MtoRref = (ffmax - ffmin) / (dpeakmax - dpeakmin);
				splineparam(MtoRref, dpeakmax, R_dcappos, dpeakmin, R_dcapneg);
				Krel = pow(MtoRref / ((Kepos + Keneg) / 2), KrR) * (Kepos + Keneg) / 2;
				Kun = fmin(KuR * MtoRref, (Kepos + Keneg) / 2);
				E = ER * (ffmax - ffmin) * ((dpeakmax - dpeakmin) - (ffmax - ffmin) / Kepos);
			}
			//------------------------------------------------------------------------------------------------------

			// -----I needed to check if my current point is above or under the new lower limit------
			if (CstateFlag != 7) {
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain); // NP
				Benchmark67_up = fspl;

				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain); // PN
				Benchmark67_down = fspl;
			}
			//---------------------------------------------------------------------------------------

			if ((ffmax_inner_inner <= Benchmark67_up && ffmax_inner_inner >= Benchmark67_down)) {

				if (CstateFlag != 7) {

					x = Cstrain;
					spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
					fouterPN = fspl;
					ekouterPN = ek;

					x = Cstrain;
					spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
					fouterNP_max = fspl;
					ekouterNP = ek;

					InCycFac = (ffmax_inner_inner - fouterPN) / (fouterNP_max - fouterPN);

				}


				x = d;
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
				fouterPN = fspl;
				ekouterPN = ek;
				x = -fabs(dpeakmax_inner_inner - d) + dpeakmax;
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
				finner = fouterPN + InCycFac * ((fouterNP_max - fabs(ffmax - fspl)) - fouterPN);
				ekinner = (ek - ekouterPN) / (((fouterNP_max - fabs(ffmax - fspl)) - fouterPN))*(InCycFac*((fouterNP_max - fabs(ffmax - fspl)) - fouterPN)) + ekouterPN;
				ek = ekinner;
				f = fmax(finner, fouterPN);
				if (f == fouterPN) {
					ek = ekouterPN;
					TstateFlag = 5;
				}
			}

			if (ffmax_inner_inner > Benchmark67_up) {

				// Before reaching the upper limit it should be a line
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d); //NP 
				fouterNP = fspl;
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d); //PN 
				fouterPN = fspl;
				ekouterPN = ek;
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmax); // PN
				f = ffmax_inner_inner - ek * (dpeakmax_inner_inner - d);
				f = fmax(f, fouterPN);
				if (f == fouterPN) {
					ek = ekouterPN;
					TstateFlag = 5;
				}
				// we should check if it has reached the upper limit or not, because after it passes the upper limit we should start the coppy-pasting process
				// also we should check that it has not gone below the lower limit

				if (f <= fouterNP && Cstress > fouterNP && f > fouterPN) {
					dpeakmax_inner_inner = d;
					ffmax_inner_inner = f;

					x = d;
					spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
					fouterPN = fspl;
					ekouterPN = ek;

					x = d;
					spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
					fouterNP_max = fspl;
					ekouterNP = ek;
					InCycFac = (ffmax_inner_inner - fouterPN) / (fouterNP_max - fouterPN);

				}
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmax); // PN
				if (f <= fouterNP && Cstress < fouterNP && f > fouterPN) {
					x = d;
					spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
					fouterPN = fspl;
					ekouterPN = ek;
					x = -fabs(dpeakmax_inner_inner - d) + dpeakmax;
					spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
					finner = fouterPN + InCycFac * ((fouterNP_max - fabs(ffmax - fspl)) - fouterPN);
					ekinner = (ek - ekouterPN) / (((fouterNP_max - fabs(ffmax - fspl)) - fouterPN))*(InCycFac*((fouterNP_max - fabs(ffmax - fspl)) - fouterPN)) + ekouterPN;
					ek = ekinner;
					f = fmax(finner, fouterPN);
					if (f == fouterPN) {
						ek = ekouterPN;
						TstateFlag = 5;
					}
				}
			}

			if (ffmax_inner_inner < Benchmark67_down) {

				if (CstateFlag != 7) {
					InCycFac = (Benchmark67_down - ffmax_inner_inner);
				}

				x = d;
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN 
				f = fspl - (InCycFac + (-InCycFac / (dpeakmax_inner_inner - dpeakmin))*(dpeakmax_inner_inner - d));
				ek = ek - (InCycFac / (dpeakmax_inner_inner - dpeakmin));
			}
			if (flagdmg == 1)
				update_damage();
			if (flagdmg_Hardening == 1)
				update_damage_hardeingin();
			if (flagdmg_Hardening_strength == 1)
				update_damage();

			checkEnvelope();
			break;

		case -4:

			define_peak();
			dirtag = 0; // Positive to Negative Direction
			MtoRref = (ffmax - ffmin) / (dpeakmax - dpeakmin);
			splineparam(MtoRref, dpeakmax, R_dcappos, dpeakmin, R_dcapneg);
			Krel = pow(MtoRref / ((Kepos + Keneg) / 2), KrR) * (Kepos + Keneg) / 2;
			Kun = fmin(KuR * MtoRref, (Kepos + Keneg) / 2);
			E = ER * (ffmax - ffmin) * ((dpeakmax - dpeakmin) - (ffmax - ffmin) / Keneg);  
			spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d);   
			f = fspl;

			if (flagdmg == 1)
				update_damage();
			if (flagdmg_Hardening == 1)
				update_damage_hardeingin();
			if (flagdmg_Hardening_strength == 1)
				update_damage();
			dpeakmax_inner = d;
			ffmax_inner = f;
			break;

		case -5:
			spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d);      //Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
			f = fspl;
			if (flagdmg == 1)
				update_damage();
			if (flagdmg_Hardening == 1)
				update_damage_hardeingin();
			if (flagdmg_Hardening_strength == 1)
				update_damage();
			checkEnvelope();
			break;

		case -6:

			define_peak();

			// Because of the damage model our peaks may be updated so we need to re-define spline parametere again 
			if (CstateFlag != -6) {
				MtoRref = (ffmax - ffmin) / (dpeakmax - dpeakmin);
				splineparam(MtoRref, dpeakmax, R_dcappos, dpeakmin, R_dcapneg);
				Krel = pow(MtoRref / ((Kepos + Keneg) / 2), KrR) * (Kepos + Keneg) / 2;
				Kun = fmin(KuR * MtoRref, (Kepos + Keneg) / 2);
				E = ER * (ffmax - ffmin) * ((dpeakmax - dpeakmin) - (ffmax - ffmin) / Keneg);
			}
			//------------------------------------------------------------------------------------------------------

			// -----I needed to check if my current point is above or under the new lower limit------
			if (CstateFlag != -6) {
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain); // NP
				Benchmark_neg_56_up = fspl;

				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain); // PN
				Benchmark_neg_56_down = fspl;
			}
			//---------------------------------------------------------------------------------------

			if ((ffmax_inner <= Benchmark_neg_56_up && ffmax_inner >= Benchmark_neg_56_down)) {

				if (CstateFlag != -6) {

					x = Cstrain;
					spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
					fouterPN = fspl;
					ekouterPN = ek;

					x = Cstrain;
					spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
					fouterNP_max = fspl;
					ekouterNP = ek;

					InCycFac = (ffmax_inner - fouterPN) / (fouterNP_max - fouterPN);

				}


				x = d;
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
				fouterPN = fspl;
				ekouterPN = ek;
				x = -fabs(dpeakmax_inner - d) + dpeakmax;
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
				finner = fouterPN + InCycFac * ((fouterNP_max - fabs(ffmax - fspl)) - fouterPN);
				ekinner = (ek - ekouterPN) / (((fouterNP_max - fabs(ffmax - fspl)) - fouterPN))*(InCycFac*((fouterNP_max - fabs(ffmax - fspl)) - fouterPN)) + ekouterPN;
				ek = ekinner;
				f = fmax(finner, fouterPN);
				if (f == fouterPN) {
					ek = ekouterPN;
					TstateFlag = 5;
				}
			}

			if (ffmax_inner > Benchmark_neg_56_up) {
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d); //NP 
				fouterNP = fspl;
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d); //PN 
				fouterPN = fspl;
				ekouterPN = ek;

				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmax); // PN
				f = ffmax_inner - ek * (dpeakmax_inner - d);
				f = fmax(f, fouterPN);
				if (f == fouterPN) {
					ek = ekouterPN;
					TstateFlag = 5;
				}
				// we should check if it has reached the upper limit or not, because after it passes the upper limit we should start the coppy-pasting process
				// also we should check that it has not gone below the lower limit
				if (f <= fouterNP && Cstress > fouterNP && f > fouterPN) {
					dpeakmax_inner = d;
					ffmax_inner = f;

					x = d;
					spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
					fouterPN = fspl;
					ekouterPN = ek;

					x = d;
					spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
					fouterNP_max = fspl;
					ekouterNP = ek;
					InCycFac = (ffmax_inner - fouterPN) / (fouterNP_max - fouterPN);

				}

				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmax); // PN

				if (f <= fouterNP && Cstress < fouterNP && f > fouterPN) {
					x = d;
					spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
					fouterPN = fspl;
					ekouterPN = ek;
					x = -fabs(dpeakmax_inner - d) + dpeakmax;
					spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
					finner = fouterPN + InCycFac * ((fouterNP_max - fabs(ffmax - fspl)) - fouterPN);
					ekinner = (ek - ekouterPN) / (((fouterNP_max - fabs(ffmax - fspl)) - fouterPN))*(InCycFac*((fouterNP_max - fabs(ffmax - fspl)) - fouterPN)) + ekouterPN;
					ek = ekinner;
					f = fmax(finner, fouterPN);
					if (f == fouterPN) {
						ek = ekouterPN;
						TstateFlag = 5;
					}
				}
			}

			if (ffmax_inner < Benchmark_neg_56_down) {

				if (CstateFlag != -6) {

					InCycFac = (Benchmark_neg_56_down - ffmax_inner);

				}
				x = d;
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN 
				f = fspl - (InCycFac + (-InCycFac / (dpeakmax_inner - dpeakmin))*(dpeakmax_inner - d));
				ek = ek - (InCycFac / (dpeakmax_inner - dpeakmin));
			}

			if (flagdmg == 1)
				update_damage();
			if (flagdmg_Hardening == 1)
				update_damage_hardeingin();
			if (flagdmg_Hardening_strength == 1)
				update_damage();

			checkEnvelope();
			break;

		case -7:
			define_peak();

			// Because of the damage model our peaks may be updated so we need to re-define spline parametere again 
			if (CstateFlag != -7) {
				MtoRref = (ffmax - ffmin) / (dpeakmax - dpeakmin);
				splineparam(MtoRref, dpeakmax, R_dcappos, dpeakmin, R_dcapneg);
				Krel = pow(MtoRref / ((Kepos + Keneg) / 2), KrR) * (Kepos + Keneg) / 2;
				Kun = fmin(KuR * MtoRref, (Kepos + Keneg) / 2);
				E = ER * (ffmax - ffmin) * ((dpeakmax - dpeakmin) - (ffmax - ffmin) / Keneg);
			}
			//------------------------------------------------------------------------------------------------------

			// -----I needed to check if my current point is above or under the new lower limit------
			if (CstateFlag != -7) {
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain); // NP
				Benchmark_neg_67_up = fspl;

				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain); // PN
				Benchmark_neg_67_down = fspl;
			}
			//---------------------------------------------------------------------------------------
			if ((ffmin_inner_inner <= Benchmark_neg_67_up && ffmin_inner_inner >= Benchmark_neg_67_down)) {

				if (CstateFlag != -7) {

					x = Cstrain;
					spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
					fouterPN_min = fspl;
					ekouterPN = ek;

					x = Cstrain;
					spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
					fouterNP = fspl;
					ekouterNP = ek;

					InCycFac = (fouterNP - ffmin_inner_inner) / (fouterNP - fouterPN_min);

				}

				x = d;
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
				fouterNP = fspl;
				ekouterNP = ek;
				x = dpeakmin + (d - dpeakmin_inner_inner);
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
				finner = fouterNP - InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)));
				ekinner = (ek - ekouterNP) / ((fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) * (InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) + ekouterNP;
				ek = ekinner;
				f = fmin(finner, fouterNP);
				if (f == fouterNP) {
					ek = ekouterNP;
					TstateFlag = -5;
				}
			}


			if (ffmin_inner_inner < Benchmark_neg_67_down) {
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d);
				fouterPN = fspl;
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d);
				fouterNP = fspl;
				ekouterNP = ek;
				spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmax); // PN
				f = ffmin_inner_inner + ek * (d - dpeakmin_inner_inner);
				f = fmin(f, fouterNP);
				if (f == fouterNP) {
					ek = ekouterNP;
					TstateFlag = -5;
				}
				// we should check if it has reached the lower limit or not, because after it passes the lower limit we should start the coppy-pasting process
				if (f >= fouterPN && Cstress < fouterPN && f < fouterNP) {
					//f = fouterNP;
					dpeakmin_inner_inner = d;
					ffmin_inner_inner = f;
					x = d;
					spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
					fouterPN_min = fspl;
					ekouterPN = ek;

					x = d;
					spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
					fouterNP = fspl;
					ekouterNP = ek;

					InCycFac = (fouterNP - ffmin_inner_inner) / (fouterNP - fouterPN_min);
				}

				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmin); // NP

				if (f >= fouterPN && Cstress > fouterPN && f < fouterNP) {
					x = d;
					spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
					fouterNP = fspl;
					ekouterNP = ek;

					x = dpeakmin + (d - dpeakmin_inner_inner);
					spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
					finner = fouterNP - InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)));
					ekinner = (ek - ekouterNP) / ((fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) * (InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) + ekouterNP;
					ek = ekinner;
					f = fmin(finner, fouterNP);
					if (f == fouterNP) {
						ek = ekouterNP;
						TstateFlag = -5;
					}
				}
			}

			if (ffmin_inner_inner > Benchmark_neg_67_up) {

				if (CstateFlag != -7) {

					InCycFac = (ffmin_inner_inner - Benchmark_neg_67_up);
				}

				x = d;
				spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP 
				f = fspl + (InCycFac + (-InCycFac / (dpeakmax - dpeakmin_inner_inner))*(d - dpeakmin_inner_inner));
				ek = ek - (InCycFac / (dpeakmax - dpeakmin_inner_inner));
			}

			if (flagdmg == 1)
				update_damage();
			if (flagdmg_Hardening == 1)
				update_damage_hardeingin();
			if (flagdmg_Hardening_strength == 1)
				update_damage();
			checkEnvelope();
			break;
		}
	}
	return 0;
}

double GMG_CyclicReinforcedConcrete::getStrain(void)
{
	return d;
}

double GMG_CyclicReinforcedConcrete::getStress(void)
{

	return f;
}

double GMG_CyclicReinforcedConcrete::getTangent(void)
{
	return ek;
}

double GMG_CyclicReinforcedConcrete::getInitialTangent(void)
{
	return (Kepos);
}

double GMG_CyclicReinforcedConcrete::getDampTangent(void)
{
	double DTangent = 0.0;
	return DTangent;
}

double GMG_CyclicReinforcedConcrete::getStrainRate(void)
{
	return 0;
}

int GMG_CyclicReinforcedConcrete::commitState(void)
{
  //commitCalledOnce = 1;

	if (TstateFlag == 5) {
		dpeakmin_inner = d;
		ffmin_inner = f;
	}

	if (TstateFlag == -5) {
		dpeakmax_inner = d;
		ffmax_inner = f;
	}

	if (TstateFlag == 6) {
		dpeakmax_inner_inner = d;
		ffmax_inner_inner = f;
	}

	if (TstateFlag == -6) {
		dpeakmin_inner_inner = d;
		ffmin_inner_inner = f;
	}

	if (TstateFlag == 7) {
		dpeakmin_inner_inner = d;
		ffmin_inner_inner = f;
	}

	if (TstateFlag == -7) {
		dpeakmax_inner_inner = d;
		ffmax_inner_inner = f;
	}

	Cstrain = d;
	Cstress = f;
	Cek = ek;
	CInCycFac = InCycFac;

	Cflagstop = flagstop;
	Cflagcap = flagcap;
	Cflagdmg = flagdmg;
	Cflagdmg_Hardening = flagdmg_Hardening;
	Cflagdmg_Hardening_strength = flagdmg_Hardening_strength;
	Ckon = kon;
	Ctrig = trig;

	Cdpeakmax = dpeakmax;
	Cdpeakmax_inner = dpeakmax_inner;
	Cdpeakmax_inner_inner = dpeakmax_inner_inner;
	Cdpealmax_bench = dpealmax_bench;
	Cdpeakmin = dpeakmin;
	Cdpeakmin_inner = dpeakmin_inner;
	Cdpeakmin_inner_inner = dpeakmin_inner_inner;
	Cdpeakmin_bench = dpeakmin_bench;
	Cffmax = ffmax;
	Cffmax_inner = ffmax_inner;
	Cffmax_inner_inner = ffmax_inner_inner;
	CfouterNP_max = fouterNP_max;
	Cffmin = ffmin;
	Cffmin_inner = ffmin_inner;
	Cffmin_inner_inner = ffmin_inner_inner;
	CfouterPN_min = fouterPN_min;

	Cspos = spos;
	Csneg = sneg;
	Cfoffpos = foffpos;
	Cfoffneg = foffneg;
	CdmgSpos = dmgSpos;
	CdmgSneg = dmgSneg;
	Calpha = alpha;
	CEt = Et;
	CEc = Ec;
	CEposnorm = Eposnorm;
	CEnegnorm = Enegnorm;

	CdeltaD = deltaD;
	CTdu = Tdu;
	CstrainRate = TstrainRate;
	CstrainMax = TstrainMax;
	CstrainMin = TstrainMin;
	CstateFlag = TstateFlag;
	Cd12 = d12;
	Cf12 = f12;
	Cd_12 = d_12;
	Cd12neg = d12neg;
	Cf12neg = f12neg;
	Cd_12neg = d_12neg;
	Cf_12neg = f_12neg;
	Cd2 = d2;
	Cf2 = f2;
	Cf_12 = f_12;
	Cd_2 = d_2;
	Cf_2 = f_2;
	Cd2neg = d2neg;
	Cf2neg = f2neg;
	Cd_2neg = d_2neg;
	Cf_2neg = f_2neg;
	Cd3 = d3;
	Cf3 = f3;
	Cd_3 = d_3;
	Cf_3 = f_3;
	Cd3neg = d3neg;
	Cf3neg = f3neg;
	Cd_3neg = d_3neg;
	Cf_3neg = f_3neg;
	Cratio = ratio;
	CdmgCounter = dmgCounter;

	CR_dcapneg = R_dcapneg;
	CR_dcappos = R_dcappos;
	Cdresneg = dresneg;
	CR_fresneg = R_fresneg;
	Cslope_neg = slope_neg;
	Cdpneg = dpneg;
	Cdrespos = drespos;
	CIntcpt_slope_neg = Intcpt_slope_neg;
	CIntcpt_slope_pos = Intcpt_slope_pos;
	CR_frespos = R_frespos;
	Cslope_pos = slope_pos;
	Cdppos = dppos;
	CR_Kdegneg = R_Kdegneg;
	CR_fcapneg = R_fcapneg;
	CR_dyneg = R_dyneg;
	CR_fyneg = R_fyneg;
	CR_Kdegpos = R_Kdegpos;
	CR_fcappos = R_fcappos;
	CR_fypos = R_fypos;
	CR_dypos = R_dypos;

	CR_dresneg = R_dresneg;
	CR_drespos = R_drespos;

	CIntcpt_deg_pos = Intcpt_deg_pos;
	CIntcpt_res_pos = Intcpt_res_pos;
	CIntcpt_deg_neg = Intcpt_deg_neg;
	CIntcpt_res_neg = Intcpt_res_neg;
	CIntcpt_Xaxis_pos = Intcpt_Xaxis_pos;
	CIntcpt_Xaxis_neg = Intcpt_Xaxis_neg;

	CInCycFac_6 = InCycFac_6;
	CInCycFac_7 = InCycFac_7;
	CInCycFac_neg_6 = InCycFac_neg_6;
	CInCycFac_neg_7 = InCycFac_neg_7;

	CBenchmark56_up = Benchmark56_up;
	CBenchmark56_down = Benchmark56_down;
	CBenchmark_neg_56_up = Benchmark_neg_56_up;
	CBenchmark_neg_56_down = Benchmark_neg_56_down;
	CBenchmark67_up = Benchmark67_up;
	CBenchmark67_down = Benchmark67_down;
	CBenchmark_neg_67_up = Benchmark_neg_67_up;
	CBenchmark_neg_67_down = Benchmark_neg_67_down;


	Cspos_pos = spos_pos;
	Csneg_pos = sneg_pos;
	CEt_pos = Et_pos;
	CEc_pos = Ec_pos;
	Cspos_neg = spos_neg;
	Csneg_neg = sneg_neg;
	CEt_neg = Et_neg;
	CEc_neg = Ec_neg;
	CEpos_pos = Epos_pos;
	CEneg_neg = Eneg_neg;

	Cspos_hard = spos_hard;
	Csneg_hard = sneg_hard;
	CEt_hard = Et_hard;
	CEc_hard = Ec_hard;
	Cspos_pos_hard = spos_pos_hard;
	Csneg_pos_hard = sneg_pos_hard;
	CEt_pos_hard = Et_pos_hard;
	CEc_pos_hard = Ec_pos_hard;
	Cspos_neg_hard = spos_neg_hard;
	Csneg_neg_hard = sneg_neg_hard;
	CEt_neg_hard = Et_neg_hard;
	CEc_neg_hard = Ec_neg_hard;
	CEpos_pos_hard = Epos_pos_hard;
	CEneg_neg_hard = Eneg_neg_hard;
	CEpos_hard = Epos_hard;
	CEneg_hard = Eneg_hard;
	CEpos = Epos;
	CEneg = Eneg;
	
	CEposnorm_hard = Eposnorm_hard;
	CEnegnorm_hard = Enegnorm_hard;
	Cdelta_pos_hard = delta_pos_hard;
	Cdelta_neg_hard = delta_neg_hard;
	Calpha_pos = alpha_pos;
	Calpha_neg = alpha_neg;

	return 0;
}

int GMG_CyclicReinforcedConcrete::revertToLastCommit(void)
{
	d = Cstrain;
	f = Cstress;
	ek = Cek;
	InCycFac = CInCycFac;

	kon = Ckon;
	trig = Ctrig;
	flagdmg = Cflagdmg;
	flagdmg_Hardening = Cflagdmg_Hardening;
	flagdmg_Hardening_strength = Cflagdmg_Hardening_strength;
	flagstop = Cflagstop;
	flagcap = Cflagcap;

	dpeakmax = Cdpeakmax;
	dpeakmax_inner = Cdpeakmax_inner;
	dpeakmax_inner_inner = Cdpeakmax_inner_inner;
	dpealmax_bench = Cdpealmax_bench;
	dpeakmin = Cdpeakmin;
	dpeakmin_inner = Cdpeakmin_inner;
	dpeakmin_inner_inner = Cdpeakmin_inner_inner;
	dpeakmin_bench = Cdpeakmin_bench;
	ffmax = Cffmax;
	ffmax_inner = Cffmax_inner;
	ffmax_inner_inner = Cffmax_inner_inner;
	fouterNP_max = CfouterNP_max;
	ffmin = Cffmin;
	ffmin_inner = Cffmin_inner;
	ffmin_inner_inner = Cffmin_inner_inner;
	fouterPN_min = CfouterPN_min;

	spos = Cspos;
	sneg = Csneg;
	foffpos = Cfoffpos;
	foffneg = Cfoffneg;
	dmgSpos = CdmgSpos;
	dmgSneg = CdmgSneg;
	alpha = Calpha;
	Et = CEt;
	Ec = CEc;
	Eposnorm = CEposnorm;
	Enegnorm = CEnegnorm;

	Tdu = CTdu;
	TstrainRate = CstrainRate;
	TstrainMax = CstrainMax;
	TstrainMin = CstrainMin;
	TstateFlag = CstateFlag;
	d12 = Cd12;
	f12 = Cf12;
	d_12 = Cd_12;
	d12neg = Cd12neg;
	f12neg = Cf12neg;
	d_12neg = Cd_12neg;
	f_12neg = Cf_12neg;
	d2 = Cd2;
	f2 = Cf2;
	f_12 = Cf_12;
	d_2 = Cd_2;
	f_2 = Cf_2;
	d2neg = Cd2neg;
	f2neg = Cf2neg;
	d_2neg = Cd_2neg;
	f_2neg = Cf_2neg;
	d3 = Cd3;
	f3 = Cf3;
	d_3 = Cd_3;
	f_3 = Cf_3;
	d3neg = Cd3neg;
	f3neg = Cf3neg;
	d_3neg = Cd_3neg;
	f_3neg = Cf_3neg;
	dmgCounter = CdmgCounter;
	R_dcapneg = CR_dcapneg;
	R_dcappos = CR_dcappos;
	dresneg = Cdresneg;
	R_fresneg = CR_fresneg;
	slope_neg = Cslope_neg;
	dpneg = Cdpneg;
	drespos = Cdrespos;
	Intcpt_slope_neg = CIntcpt_slope_neg;
	Intcpt_slope_pos = CIntcpt_slope_pos;
	R_frespos = CR_frespos;
	slope_pos = Cslope_pos;
	dppos = Cdppos;
	R_Kdegneg = CR_Kdegneg;
	R_fcapneg = CR_fcapneg;
	R_dyneg = CR_dyneg;
	R_fyneg = CR_fyneg;
	R_Kdegpos = CR_Kdegpos;
	R_fcappos = CR_fcappos;
	R_dypos = CR_dypos;
	R_fypos = CR_fypos;

	R_dresneg = CR_dresneg;
	R_drespos = CR_drespos;
	Intcpt_deg_pos = CIntcpt_deg_pos;
	Intcpt_res_pos = CIntcpt_res_pos;
	Intcpt_deg_neg = CIntcpt_deg_neg;
	Intcpt_res_neg = CIntcpt_res_neg;
	Intcpt_Xaxis_pos = CIntcpt_Xaxis_pos;
	Intcpt_Xaxis_neg = CIntcpt_Xaxis_neg;

	InCycFac_6 = CInCycFac_6;
	InCycFac_7 = CInCycFac_7;
	InCycFac_neg_6 = CInCycFac_neg_6;
	InCycFac_neg_7 = CInCycFac_neg_7;

	Benchmark56_down = CBenchmark56_down;
	Benchmark56_up = CBenchmark56_up;
	Benchmark_neg_56_up = CBenchmark_neg_56_up;
	Benchmark_neg_56_down = CBenchmark_neg_56_down;
	Benchmark67_up = CBenchmark67_up;
	Benchmark67_down = CBenchmark67_down;
	Benchmark_neg_67_up = CBenchmark_neg_67_up;
	Benchmark_neg_67_down = CBenchmark_neg_67_down;


	spos_pos = Cspos_pos;
	sneg_pos = Csneg_pos;
	Et_pos = CEt_pos;
	Ec_pos = CEc_pos;
	spos_neg = Cspos_neg;
	sneg_neg = Csneg_neg;
	Et_neg = CEt_neg;
	Ec_neg = CEc_neg;
	Epos_pos = CEpos_pos;
	Eneg_neg = CEneg_neg;
	T_area = CT_area;


	spos_hard = Cspos_hard;
	sneg_hard = Csneg_hard;
	Et_hard = CEt_hard;
	Ec_hard = CEc_hard;
	spos_pos_hard = Cspos_pos_hard;
	sneg_pos_hard = Csneg_pos_hard;
	Et_pos_hard = CEt_pos_hard;
	Ec_pos_hard = CEc_pos_hard;
	spos_neg_hard = Cspos_neg_hard;
	sneg_neg_hard = Csneg_neg_hard;
	Et_neg_hard = CEt_neg_hard;
	Ec_neg_hard = CEc_neg_hard;
	Epos_pos_hard = CEpos_pos_hard;
	Eneg_neg_hard = CEneg_neg_hard;
	Epos_hard = CEpos_hard;
	Eneg_hard = CEneg_hard;
	Epos = CEpos;
	Eneg = CEneg;
	
	Eposnorm_hard = CEposnorm_hard;
	Enegnorm_hard = CEnegnorm_hard;
	delta_pos_hard = Cdelta_pos_hard;
	delta_neg_hard = Cdelta_neg_hard;
	alpha_pos = Calpha_pos;
	alpha_neg = Calpha_neg;

	return 0;
}

int GMG_CyclicReinforcedConcrete::revertToStart(void)
{
  //commitCalledOnce = 0;
	d = Cstrain = 0.0;
	f = Cstress = 0.0;
	ek = Cek = Kepos;

	kon = Ckon = 0;
	trig = Ctrig = 0;
	flagstop = Cflagstop = 0;
	flagcap = Cflagcap = 0;

	dypos = fypos / Kepos;
	drespos = dcappos + (frespos - fcappos) / Kdegpos;
	dUpos = drespos + delUpos;
	dppos = dcappos - dypos;
	dpcpos = drespos - dcappos;
	ahardpos = ((fcappos - fypos) / dppos) / Kepos;
	adegpos = Kdegpos / Kepos;

	dyneg = fyneg / Keneg;
	dresneg = dcapneg + (fresneg - fcapneg) / Kdegneg;
	dUneg = dresneg + delUneg;
	dpneg = dcapneg - dyneg;
	dpcneg = dresneg - dcapneg;
	ahardneg = ((fcapneg - fyneg) / dpneg) / Keneg;
	adegneg = Kdegneg / Keneg;

	dpeakmax = Cdpeakmax = dpealmax_bench = Cdpealmax_bench = dypos;
	dpeakmin = Cdpeakmin = dpeakmin_bench = Cdpeakmin_bench = dyneg;
	ffmax = Cffmax = fypos;
	ffmin = Cffmin = fyneg;

	foffpos = fcappos - Kdegpos * dcappos;
	Ed0pos = 0.0;
	Ed1pos = Ed0pos + fabs(0.5 * fypos * dypos) + (0.5 * fabs(fypos + fcappos) * fabs(dppos)) + (0.5 * fabs(fcappos + frespos) * fabs(dpcpos)) + 0.5 * fabs(frespos) * fabs(frespos / Kdegpos);
	Ed0pos_hard = 0.0;
	Ed1pos_hard = Ed0pos_hard + abs(fypos * dypos * (Ref_Energy_Coe - 1)) + abs(fypos * dypos) / 2.0;
	dmgSpos = CdmgSpos = 0.0;
	spos = Cspos = 0.0;
	Et = CEt = 0.0;
	Eposnorm = CEposnorm = 0.0;

	foffneg = fcapneg - Kdegneg * dcapneg;
	Ed0neg = 0.0;
	Ed1neg = Ed0neg + (0.5 * fabs(fyneg * dyneg)) + (0.5 * fabs(fyneg + fcapneg) * fabs(dpneg)) + (0.5 * fabs(fcapneg + fresneg) * fabs(dpcneg)) + 0.5 * fabs(fresneg) * fabs(fresneg / Kdegneg);
	Ed0neg_hard = 0.0;
	Ed1neg_hard = Ed0neg_hard + (abs(fyneg * dyneg * (Ref_Energy_Coe - 1)) + abs(fyneg * dyneg) / 2.0);
	dmgSneg = CdmgSneg = 0.0;
	alpha = Calpha = 0.0;
	sneg = Csneg = 0.0;
	Ec = CEc = 0.0;
	Enegnorm = CEnegnorm = 0.0;

	//---------------------------My Parameters-----------------------
	CstateFlag = 0;
	flagdmg = 0;
	flagdmg_Hardening = 0;
	flagdmg_Hardening_strength = Cflagdmg_Hardening_strength = 0;
	Et = 0.0;
	Ec = 0.0;
	dmgCounter = CdmgCounter = 0;
	spos_pos = Cspos_pos = 0.0;
	sneg_pos = Csneg_pos = 0.0;
	Et_pos = CEt_pos = 0.0;
	Ec_pos = CEc_pos = 0.0;
	spos_neg = Cspos_neg = 0.0;
	sneg_neg = Csneg_neg = 0.0;
	Et_neg = CEt_neg = 0.0;
	Ec_neg = CEc_neg = 0.0;
	Epos_pos = CEpos_pos = 0.0;
	Eneg_neg = CEneg_neg = 0.0;
	T_area = CT_area = 0.0;

	spos_hard = Cspos_hard = 0.0;
	sneg_hard = Csneg_hard = 0.0;
	Et_hard = CEt_hard = 0.0;
	Ec_hard = CEc_hard = 0.0;
	spos_pos_hard = Cspos_pos_hard = 0.0;
	sneg_pos_hard = Csneg_pos_hard = 0.0;
	Et_pos_hard = CEt_pos_hard = 0.0;
	Ec_pos_hard = CEc_pos_hard = 0.0;
	spos_neg_hard = Cspos_neg_hard = 0.0;
	sneg_neg_hard = Csneg_neg_hard = 0.0;
	Et_neg_hard = CEt_neg_hard = 0.0;
	Ec_neg_hard = CEc_neg_hard = 0.0;
	Epos_pos_hard = CEpos_pos_hard = 0.0;
	Eneg_neg_hard = CEneg_neg_hard = 0.0;
	CEpos_hard = Epos_hard = 0.0;
	CEneg_hard = Eneg_hard = 0.0;
	CEpos = Epos = 0.0;
	CEneg = Eneg = 0.0;
	CEposnorm_hard = Eposnorm_hard = 0.0;
	CEnegnorm_hard = Enegnorm_hard = 0.0;
	Cdelta_pos_hard = delta_pos_hard = 0.0;
	Cdelta_neg_hard = delta_neg_hard = 0.0;
	Calpha_pos = alpha_pos = 0.0;
	Calpha_neg = alpha_neg = 0.0;
	return 0;
}

UniaxialMaterial *
GMG_CyclicReinforcedConcrete::getCopy(void)
{
	GMG_CyclicReinforcedConcrete *theCopy = new GMG_CyclicReinforcedConcrete(this->getTag(), Kepos, Keneg, fypos, fyneg, fcappos, fcapneg,
		dcappos, dcapneg, Kdegpos, Kdegneg,
		frespos, fresneg, delUpos, delUneg,
		alpha_Er_Post_Yielding, beta_Er_Post_Yielding,
		alpha_Er_Post_Capping, beta_Er_Post_Capping,
		Er_Post_Yielding, Er_Post_Capping,
		Kun_Post_Yielding, Kun_Post_Capping, Kr_Post_Yielding, Kr_Post_Capping,
		delta_ratio_max_hard, Ref_Energy_Coe,
		C1, C2, C3, solpe_post_yielding, solpe_post_capping);

	*theCopy = *this;
	
	return theCopy;
}

int
GMG_CyclicReinforcedConcrete::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	static Vector data(180);
	data(0) = this->getTag();

	// Monotonic Backbone
	data(1) = Kepos;
	data(2) = Keneg;
	data(3) = fypos;
	data(4) = fyneg;
	data(5) = fcappos;
	data(6) = fcapneg;
	data(7) = dcappos;
	data(8) = dcapneg;
	data(9) = Kdegpos;
	data(10) = Kdegneg;
	data(11) = frespos;
	data(12) = fresneg;
	data(13) = delUpos;
	data(14) = delUneg;

	// Cyclic Properties
	data(15) = alpha_Er_Post_Yielding;
	data(16) = beta_Er_Post_Yielding;
	data(17) = alpha_Er_Post_Capping;
	data(18) = beta_Er_Post_Capping;
	data(19) = Er_Post_Yielding;
	data(20) = Er_Post_Capping;
	data(21) = Kun_Post_Yielding;
	data(22) = Kun_Post_Capping;
	data(23) = Kr_Post_Yielding;
	data(24) = Kr_Post_Capping;

	// Damage Properties
	data(25) = C1;
	data(26) = C2;
	data(27) = C3;
	data(28) = solpe_post_capping;
	data(29) = delta_ratio_max_hard;
	data(30) = Ref_Energy_Coe;
	data(31) = solpe_post_yielding;

	// Converged Variables
	data(32) = Cstrain;
	data(33) = Cstress;
	data(34) = Cek;
	data(35) = CInCycFac;

	data(36) = Cflagstop;
	data(37) = Cflagcap;
	data(38) = Cflagdmg;
	data(39) = Cflagdmg_Hardening;
	data(40) = Cflagdmg_Hardening_strength;
	data(41) = Ckon;
	data(42) = Ctrig;

	data(43) = Cdpeakmax;
	data(44) = Cdpeakmax_inner;
	data(45) = Cdpeakmax_inner_inner;
	data(46) = Cdpealmax_bench;
	data(47) = Cdpeakmin;
	data(48) = Cdpeakmin_inner;
	data(49) = Cdpeakmin_inner_inner;
	data(50) = Cdpeakmin_bench;
	data(51) = Cffmax;
	data(52) = Cffmax_inner;
	data(53) = Cffmax_inner_inner;
	data(54) = CfouterNP_max;
	data(55) = Cffmin;
	data(56) = Cffmin_inner;
	data(57) = Cffmin_inner_inner;
	data(58) = CfouterPN_min;

	data(59) = Cspos;
	data(60) = Csneg;
	data(61) = Cfoffpos;
	data(62) = Cfoffneg;
	data(63) = CdmgSpos;
	data(64) = CdmgSneg;
	data(65) = Calpha;	
	data(66) = CEt;
	data(67) = CEc;
	data(68) = CEposnorm;
	data(69) = CEnegnorm;

	data(70) = CdeltaD;
	data(71) = CTdu;
	data(72) = CstrainRate;
	data(73) = CstrainMax;
	data(74) = CstrainMin;
	data(75) = CstateFlag;
	data(76) = Cd12;
	data(77) = Cf12;
	data(78) = Cd_12;
	data(79) = Cf_12;
	data(80) = Cd12neg;
	data(81) = Cf12neg;
	data(82) = Cd_12neg;
	data(83) = Cf_12neg;	
	data(84) = Cd2;
	data(85) = Cf2;
	data(86) = Cd_2;
	data(87) = Cf_2;	
	data(88) = Cd2neg;
	data(89) = Cf2neg;
	data(90) = Cd_2neg;
	data(91) = Cf_2neg;
	data(92) = Cd3;
	data(93) = Cf3;
	data(94) = Cd_3;
	data(95) = Cf_3;	
	data(96) = Cd3neg;
	data(97) = Cf3neg;
	data(98) = Cd_3neg;
	data(99) = Cf_3neg;	
	data(100) = CdmgCounter;

	data(101) = CR_dcapneg;
	data(102) = CR_dcappos;	
	data(103) = Cdresneg;
	data(104) = CR_fresneg;
	data(105) = Cslope_neg;
	data(106) = Cdpneg;
	data(107) = Cdrespos;
	data(108) = CIntcpt_slope_neg;
	data(109) = CIntcpt_slope_pos;
	data(110) = CR_frespos;
	data(111) = Cslope_pos;
	data(112) = Cdppos;
	data(113) = CR_Kdegneg;
	data(114) = CR_fcapneg;
	data(115) = CR_dyneg;
	data(116) = CR_fyneg;
	data(117) = CR_Kdegpos;
	data(118) = CR_fcappos;
	data(119) = CR_fypos;
	data(120) = CR_dypos;	
	data(121) = CR_dresneg;
	data(122) = CR_drespos;

	data(123) = CIntcpt_deg_pos;
	data(124) = CIntcpt_res_pos;
	data(125) = CIntcpt_deg_neg;	
	data(126) = CIntcpt_res_neg;
	data(127) = CIntcpt_Xaxis_pos;
	data(128) = CIntcpt_Xaxis_neg;
	
	data(129) = CInCycFac_6;
	data(130) = CInCycFac_7;
	data(131) = CInCycFac_neg_6;
	data(132) = CInCycFac_neg_7;

	data(133) = CBenchmark56_up;
	data(134) = CBenchmark56_down;
	data(135) = CBenchmark_neg_56_up;
	data(136) = CBenchmark_neg_56_down;
	data(137) = CBenchmark67_up;
	data(138) = CBenchmark67_down;
	data(139) = CBenchmark_neg_67_up;
	data(140) = CBenchmark_neg_67_down;

	data(141) = Cspos_pos;
	data(142) = Csneg_pos;
	data(143) = CEt_pos;
	data(144) = CEc_pos;
	data(145) = Cspos_neg;
	data(146) = Csneg_neg;
	data(147) = CEt_neg;
	data(148) = CEc_neg;
	data(149) = CEpos_pos;
	data(150) = CEneg_neg;
	data(151) = CT_area;

	data(152) = Cspos_hard;
	data(153) = Csneg_hard;
	data(154) = CEt_hard;
	data(155) = CEc_hard;
	data(156) = Cspos_pos_hard;
	data(157) = Csneg_pos_hard;
	data(158) = CEt_pos_hard;
	data(159) = CEc_pos_hard;
	data(160) = Cspos_neg_hard;
	data(161) = Csneg_neg_hard;
	data(162) = CEt_neg_hard;
	data(163) = CEc_neg_hard;
	data(164) = CEpos_pos_hard;
	data(165) = CEneg_neg_hard;
	data(166) = CEpos_hard;
	data(167) = CEneg_hard;
	data(168) = CEpos;
	data(169) = CEneg;
	
	data(170) = CEposnorm_hard;
	data(171) = CEnegnorm_hard;
	data(172) = Cdelta_pos_hard;
	data(173) = Cdelta_neg_hard;
	data(174) = Calpha_pos;
	data(175) = Calpha_neg;

	data(176) = Ed0pos_hard;
	data(177) = Ed1pos_hard;
	data(178) = Ed0neg_hard;
	data(179) = Ed1neg_hard;

	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "GMG_CyclicReinforcedConcrete::sendSelf() - failed to send data\n";
	return res;
}

int
GMG_CyclicReinforcedConcrete::recvSelf(int cTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;
	static Vector data(180);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);

	if (res < 0) {
		opserr << "GMG_CyclicReinforcedConcrete::recvSelf() - failed to receive data\n";
		this->setTag(0);
	}
	else {
		this->setTag((int)data(0));
		// Monotonic Backbone
		Kepos = data(1);
		Keneg = data(2);
		fypos = data(3);
		fyneg = data(4);
		fcappos = data(5);
		fcapneg = data(6);
		dcappos = data(7);
		dcapneg = data(8);
		Kdegpos = data(9);
		Kdegneg = data(10);
		frespos = data(11);
		fresneg = data(12);
		delUpos = data(13);
		delUneg = data(14);

		// Cyclic Properties
		alpha_Er_Post_Yielding = data(15);
		beta_Er_Post_Yielding = data(16);
		alpha_Er_Post_Capping = data(17);
		beta_Er_Post_Capping = data(18);
		Er_Post_Yielding = data(19);
		Er_Post_Capping = data(20);
		Kun_Post_Yielding = data(21);
		Kun_Post_Capping = data(22);
		Kr_Post_Yielding = data(23);
		Kr_Post_Capping = data(24);

		// Damage Properties
		C1 = data(25);
		C2 = data(26);
		C3 = data(27);
		solpe_post_capping = data(28);
		delta_ratio_max_hard = data(29);
		Ref_Energy_Coe = data(30);
		solpe_post_yielding = data(31);

		// Converged Variables
		Cstrain = data(32);
		Cstress = data(33);
		Cek = data(34);
		CInCycFac = data(35);

		Cflagstop = data(36);
		Cflagcap = data(37);
		Cflagdmg = data(38);
		Cflagdmg_Hardening = data(39);
		Cflagdmg_Hardening_strength = data(40);
		Ckon = data(41);
		Ctrig = data(42);

		Cdpeakmax = data(43);
		Cdpeakmax_inner = data(44);
		Cdpeakmax_inner_inner = data(45);
		Cdpealmax_bench = data(46);
		Cdpeakmin = data(47);
		Cdpeakmin_inner = data(48);
		Cdpeakmin_inner_inner = data(49);
		Cdpeakmin_bench = data(50);
		Cffmax = data(51);
		Cffmax_inner = data(52);
		Cffmax_inner_inner = data(53);
		CfouterNP_max = data(54);
		Cffmin = data(55);
		Cffmin_inner = data(56);
		Cffmin_inner_inner = data(57);
		CfouterPN_min = data(58);

		Cspos = data(59);
		Csneg = data(60);
		Cfoffpos = data(61);
		Cfoffneg = data(62);
		CdmgSpos = data(63);
		CdmgSneg = data(64);
		Calpha = data(65);	
		CEt = data(66);
		CEc = data(67);
		CEposnorm = data(68);
		CEnegnorm = data(69);

		CdeltaD = data(70);
		CTdu = data(71);
		CstrainRate = data(72);
		CstrainMax = data(73);
		CstrainMin = data(74);
		CstateFlag = data(75);
		Cd12 = data(76);
		Cf12 = data(77);
		Cd_12 = data(78);
		Cf_12 = data(79);
		Cd12neg = data(80);
		Cf12neg = data(81);
		Cd_12neg = data(82);
		Cf_12neg = data(83);
		Cd2 = data(84);
		Cf2 = data(85);
		Cd_2 = data(86);
		Cf_2 = data(87);
		Cd2neg = data(88);
		Cf2neg = data(89);
		Cd_2neg = data(90);
		Cf_2neg = data(91);
		Cd3 = data(92);
		Cf3 = data(93);
		Cd_3 = data(94);
		Cf_3 = data(95);
		Cd3neg = data(96);
		Cf3neg = data(97);
		Cd_3neg = data(98);
		Cf_3neg = data(99);
		CdmgCounter = data(100);
		
		CR_dcapneg = data(101);
		CR_dcappos = data(102);
		Cdresneg = data(103);
		CR_fresneg = data(104);
		Cslope_neg = data(105);
		Cdpneg = data(106);
		Cdrespos = data(107);
		CIntcpt_slope_neg = data(108);
		CIntcpt_slope_pos = data(109);
		CR_frespos = data(110);
		Cslope_pos = data(111);
		Cdppos = data(112);
		CR_Kdegneg = data(113);
		CR_fcapneg = data(114);
		CR_dyneg = data(115);
		CR_fyneg = data(116);
		CR_Kdegpos = data(117);
		CR_fcappos = data(118);
		CR_fypos = data(119);
		CR_dypos = data(120);
		CR_dresneg = data(121);
		CR_drespos = data(122);
		
		CIntcpt_deg_pos = data(123);
		CIntcpt_res_pos = data(124);
		CIntcpt_deg_neg = data(125);	
		CIntcpt_res_neg = data(126);
		CIntcpt_Xaxis_pos = data(127);
		CIntcpt_Xaxis_neg = data(128);
		
		CInCycFac_6 = data(129);
		CInCycFac_7 = data(130);
		CInCycFac_neg_6 = data(131);
		CInCycFac_neg_7 = data(132);
		
		CBenchmark56_up = data(133);
		CBenchmark56_down = data(134);
		CBenchmark_neg_56_up = data(135);
		CBenchmark_neg_56_down = data(136);
		CBenchmark67_up = data(137);
		CBenchmark67_down = data(138);
		CBenchmark_neg_67_up = data(139);
		CBenchmark_neg_67_down = data(140);
		
		Cspos_pos = data(141);
		Csneg_pos = data(142);
		CEt_pos = data(143);
		CEc_pos = data(144);
		Cspos_neg = data(145);
		Csneg_neg = data(146);
		CEt_neg = data(147);
		CEc_neg = data(148);
		CEpos_pos = data(149);
		CEneg_neg = data(150);
		CT_area = data(151);
		
		Cspos_hard = data(152);
		Csneg_hard = data(153);
		CEt_hard = data(154);
		CEc_hard = data(155);
		Cspos_pos_hard = data(156);
		Csneg_pos_hard = data(157);
		CEt_pos_hard = data(158);
		CEc_pos_hard = data(159);
		Cspos_neg_hard = data(160);
		Csneg_neg_hard = data(161);
		CEt_neg_hard = data(162);
		CEc_neg_hard = data(163);
		CEpos_pos_hard = data(164);
		CEneg_neg_hard = data(165);
		CEpos_hard = data(166);
		CEneg_hard = data(167);
		CEpos = data(168);
		CEneg = data(169);
		
		CEposnorm_hard = data(170);
		CEnegnorm_hard = data(171);
		Cdelta_pos_hard = data(172);
		Cdelta_neg_hard = data(173);
		Calpha_pos = data(174);
		Calpha_neg = data(175);
		
		Ed0pos_hard = data(176);
		Ed1pos_hard = data(177);
		Ed0neg_hard = data(178);
		Ed1neg_hard = data(179);
		  
		this->revertToLastCommit();
	}

	return res;
}

void
GMG_CyclicReinforcedConcrete::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "GMG_CyclicReinforcedConcrete tag: " << this->getTag() << endln;
		s << "  Kepos: " << Kepos << endln;
		s << "  Keneg: " << Keneg << endln;
		s << "  fypos: " << fypos << endln;
		s << "  fyneg: " << fyneg << endln;
		s << "  fcappos: " << fcappos << endln;
		s << "  fcapneg: " << fcapneg << endln;
		s << "  dcappos: " << dcappos << endln;
		s << "  dcapneg: " << dcapneg << endln;
		s << "  Kdegpos: " << Kdegpos << endln;
		s << "  Kdegneg: " << Kdegneg << endln;
		s << "  frespos: " << frespos << endln;
		s << "  fresneg: " << fresneg << endln;
		s << "  delUpos: " << delUpos << endln;
		s << "  delUneg: " << delUneg << endln;
		s << "  alpha_Er_Post_Yielding: " << alpha_Er_Post_Yielding << endln;
		s << "  beta_Er_Post_Yielding: " << beta_Er_Post_Yielding << endln;
		s << "  alpha_Er_Post_Capping: " << alpha_Er_Post_Capping << endln;
		s << "  beta_Er_Post_Capping: " << beta_Er_Post_Capping << endln;
		s << "  Er_Post_Yielding: " << Er_Post_Yielding << endln;
		s << "  Er_Post_Capping: " << Er_Post_Capping << endln;
		s << "  Kun_Post_Yielding: " << Kun_Post_Yielding << endln;
		s << "  Kun_Post_Capping: " << Kun_Post_Capping << endln;
		s << "  Kr_Post_Yielding: " << Kr_Post_Yielding << endln;
		s << "  Kr_Post_Capping: " << Kr_Post_Capping << endln;
		s << "  delta_ratio_max_hard: " << delta_ratio_max_hard << endln;
		s << "  Ref_Energy_Coe: " << Ref_Energy_Coe << endln;
		s << "  C1: " << C1 << endln;
		s << "  C2: " << C2 << endln;
		s << "  C3: " << C3 << endln;
		s << "  solpe_post_yielding: " << solpe_post_yielding << endln;
		s << "  solpe_post_capping: " << solpe_post_capping << endln;
	}



	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"GMG_CyclicReinforcedConcrete\", ";
		s << "\"Kepos\": " << Kepos << ", ";
		s << "\"Keneg\": " << Keneg << ", ";
		s << "\"fypos\": " << fypos << ", ";
		s << "\"fyneg\": " << fyneg << ", ";
		s << "\"fcappos\": " << fcappos << ", ";
		s << "\"fcapneg\": " << fcapneg << ", ";
		s << "\"dcappos\": " << dcappos << ", ";
		s << "\"dcapneg\": " << dcapneg << ", ";
		s << "\"Kdegpos\": " << Kdegpos << ", ";
		s << "\"Kdegneg\": " << Kdegneg << ", ";
		s << "\"frespos\": " << frespos << ", ";
		s << "\"fresneg\": " << fresneg << ", ";
		s << "\"delUpos\": " << delUpos << ", ";
		s << "\"delUneg\": " << delUneg << ", ";
		s << "\"alpha_Er_Post_Yielding\": " << alpha_Er_Post_Yielding << ", ";
		s << "\"beta_Er_Post_Yielding\": " << beta_Er_Post_Yielding << ", ";
		s << "\"alpha_Er_Post_Capping\": " << alpha_Er_Post_Capping << ", ";
		s << "\"beta_Er_Post_Capping\": " << beta_Er_Post_Capping << ", ";
		s << "\"Er_Post_Yielding\": " << Er_Post_Yielding << ", ";
		s << "\"Er_Post_Capping\": " << Er_Post_Capping << ", ";
		s << "\"Kun_Post_Yielding\": " << Kun_Post_Yielding << ", ";
		s << "\"Kun_Post_Capping\": " << Kun_Post_Capping << ", ";
		s << "\"Kr_Post_Yielding\": " << Kr_Post_Yielding << ", ";
		s << "\"Kr_Post_Capping\": " << Kr_Post_Capping << ", ";
		s << "\"delta_ratio_max_hard\": " << delta_ratio_max_hard << ", ";
		s << "\"Ref_Energy_Coe\": " << Ref_Energy_Coe << ", ";
		s << "\"C1\": " << C1 << ", ";
		s << "\"C2\": " << C2 << ", ";
		s << "\"C3\": " << C3 << ", ";
		s << "\"solpe_post_yielding\": " << solpe_post_yielding << ", ";
		s << "\"solpe_post_capping\": " << solpe_post_capping << ", ";
	}
}


void
GMG_CyclicReinforcedConcrete::defineBackbone(void)
{

	BenMark = d;
	if (BenMark >= 0) {
		R_dypos = fypos / Kepos;//ok ok
		R_fypos = fypos;//ok ok
		R_dcappos = dcappos;//ok ok
		R_fcappos = fcappos;//ok ok
		R_Kdegpos = Kdegpos; //ok ok

		R_fyneg = fyneg; //ok ok
		R_dyneg = fyneg / Keneg; //ok ok
		R_fcapneg = fcapneg;//ok ok
		R_dcapneg = dcapneg;//ok ok
		R_Kdegneg = Kdegneg;//ok ok

		dppos = dcappos - dypos; //ok ok
		slope_pos = (fcappos - fypos) / dppos; //ok ok
		R_frespos = frespos; //ok ok
		drespos = dcappos + (frespos - fcappos) / Kdegpos; //ok ok
		R_drespos = drespos; //ok ok
		Intcpt_slope_pos = fcappos - (slope_pos * R_dcappos); //ok ok
		Intcpt_deg_pos = fabs(R_fcappos - Kdegpos * R_dcappos); //ok ok
		Intcpt_res_pos = fabs(R_frespos - Kdegpos * delUpos); //ok ok
		Intcpt_Xaxis_pos = (Kdegpos * delUpos - R_frespos) / Kdegpos; //ok ok

		dpneg = R_dcapneg - R_dyneg; //ok ok
		slope_neg = (R_fcapneg - R_fyneg) / dpneg; //ok ok
		R_fresneg = fresneg; //ok ok
		dresneg = R_dcapneg + (R_fresneg - R_fcapneg) / R_Kdegneg; //ok ok
		R_dresneg = dresneg; //ok ok
		Intcpt_slope_neg = R_fcapneg - (slope_neg * R_dcapneg); //ok ok
		Intcpt_deg_neg = fabs(R_fcapneg - Kdegneg * R_dcapneg); //ok ok
		Intcpt_res_neg = fabs(R_fresneg - Kdegneg * delUneg); //ok ok
		Intcpt_Xaxis_neg = (Kdegneg * delUneg - R_fresneg) / Kdegneg; //ok ok
	}

	else if (BenMark < 0) {
		R_dyneg = fyneg / Keneg;
		R_fyneg = fyneg;
		R_dcapneg = dcapneg;
		R_fcapneg = fcapneg;
		R_Kdegneg = Kdegneg;

		R_fypos = fypos;
		R_dypos = fypos / Kepos;
		R_dcappos = dcappos;
		R_fcappos = fcappos;

		dpneg = R_dcapneg - R_dyneg;
		slope_neg = (R_fcapneg - R_fyneg) / dpneg;
		R_fresneg = fresneg;
		dresneg = R_dcapneg + (R_fresneg - R_fcapneg) / R_Kdegneg;
		R_dresneg = dresneg;
		Intcpt_slope_neg = R_fcapneg - (slope_neg * R_dcapneg);
		Intcpt_deg_neg = fabs(R_fcapneg - Kdegneg * R_dcapneg);
		Intcpt_res_neg = fabs(R_fresneg - Kdegneg * delUneg);
		Intcpt_Xaxis_neg = (Kdegneg * delUneg - R_fresneg) / Kdegneg;

		dppos = R_dcappos - R_dypos;
		slope_pos = (R_fcappos - R_fypos) / dppos;
		R_frespos = frespos;
		drespos = dcappos + (frespos - fcappos) / Kdegpos;
		R_drespos = drespos;
		Intcpt_slope_pos = R_fcappos - (slope_pos * R_dcappos);
		Intcpt_deg_pos = fabs(R_fcappos - Kdegpos * R_dcappos);
		Intcpt_res_pos = fabs(R_frespos - Kdegpos * delUpos);
		Intcpt_Xaxis_pos = (Kdegpos * delUpos - R_frespos) / Kdegpos;
	}

}

int
GMG_CyclicReinforcedConcrete::getStateFlag(void)
{

	if ((CstateFlag == 0 || CstateFlag == 12 || CstateFlag == -5 || CstateFlag == 6 || CstateFlag == -7 || CstateFlag == 5 || CstateFlag == 4) && Tdu > 0.0 /*&& d >= dypos*/ && d >= dpeakmax && d < R_dcappos && BenMark > 0)
		return 12;
	else if ((CstateFlag == 0 || CstateFlag == 12 || CstateFlag == -5 || CstateFlag == 6 || CstateFlag == -7 || CstateFlag == 5 || CstateFlag == 4) && Tdu > 0.0 /*&& d >= dypos*/ && d >= dpeakmax && d < R_dcappos && BenMark < 0)
		return 12;
	else if ((CstateFlag == 12 || CstateFlag == 2 || CstateFlag == -5 || CstateFlag == 6 || CstateFlag == -7 || CstateFlag == 5 || CstateFlag == 4) && Tdu > 0.0 && d >= R_dcappos && d >= dpeakmax && d < R_drespos && BenMark > 0)
		return 2;
	else if ((CstateFlag == 12 || CstateFlag == 2 || CstateFlag == -5 || CstateFlag == 6 || CstateFlag == -7 || CstateFlag == 5 || CstateFlag == 4) && Tdu > 0.0 && d >= R_dcappos && d >= dpeakmax && d < R_drespos && BenMark < 0)
		return 2;
	else if ((CstateFlag == 12 || CstateFlag == 2 || CstateFlag == 3 || CstateFlag == -5 || CstateFlag == 6 || CstateFlag == -7 || CstateFlag == 5 || CstateFlag == 4) && Tdu > 0.0 && d >= R_drespos && d < delUpos && d >= dpeakmax && BenMark > 0)
		return 3;
	else if ((CstateFlag == 12 || CstateFlag == 2 || CstateFlag == 3 || CstateFlag == -5 || CstateFlag == 6 || CstateFlag == -7 || CstateFlag == 5 || CstateFlag == 4) && Tdu > 0.0 && d >= R_drespos && d < delUpos && d >= dpeakmax && BenMark < 0)
		return 3;

	else if ((CstateFlag == 3 || CstateFlag == 30 || CstateFlag == 31) && Tdu > 0.0 && d >= delUpos && d <= Intcpt_Xaxis_pos && d >= dpeakmax && BenMark > 0)
		return 30;
	else if ((CstateFlag == 3 || CstateFlag == 30 || CstateFlag == 31) && Tdu > 0.0 && d >= delUpos && d <= Intcpt_Xaxis_pos && d >= dpeakmax && BenMark < 0)
		return 30;
	else if ((CstateFlag == 30 || CstateFlag == 40) && Tdu > 0.0 && d >= Intcpt_Xaxis_pos && d >= dpeakmax && BenMark > 0)
		return 40;
	else if ((CstateFlag == 30 || CstateFlag == 40) && Tdu > 0.0 && d >= Intcpt_Xaxis_pos && d >= dpeakmax && BenMark < 0)
		return 40;
	else if ((CstateFlag == 30 || CstateFlag == 31) && (Tdu > 0.0 || Tdu <= 0.0) && (d > dpeakmin || d < dpeakmax))
		return 31;
	else if ((CstateFlag == 40 || CstateFlag == 41) && (Tdu > 0.0 || Tdu <= 0.0))
		return 41;

	else if ((CstateFlag == 12 || CstateFlag == 2 || CstateFlag == 3) && Tdu < 0.0)
		return 4;
	else if ((CstateFlag == 4 || CstateFlag == 5) && Tdu < 0.0 && d > dpeakmin)
		return 5;

	else if ((CstateFlag == 4 || CstateFlag == 5 || CstateFlag == 6) && Tdu > 0.0 && d < dpeakmax)
		return 6;
	else if ((CstateFlag == 6 || CstateFlag == 7 || CstateFlag == -7) && Tdu < 0.0 && d > dpeakmin)
		return 7;

	else if ((CstateFlag == 0 || CstateFlag == -12 || CstateFlag == 5 || CstateFlag == -6 || CstateFlag == 7 || CstateFlag == -5 || CstateFlag == -4) && Tdu < 0.0 && d <= dpeakmin && d > R_dcapneg && BenMark > 0)
		return -12;
	else if ((CstateFlag == 0 || CstateFlag == -12 || CstateFlag == 5 || CstateFlag == -6 || CstateFlag == 7 || CstateFlag == -5 || CstateFlag == -4) && Tdu < 0.0 && d <= dpeakmin && d > R_dcapneg && BenMark < 0)
		return -12;
	else if ((CstateFlag == -12 || CstateFlag == -2 || CstateFlag == 5 || CstateFlag == -6 || CstateFlag == 7 || CstateFlag == -5 || CstateFlag == -4) && Tdu < 0.0 && d <= R_dcapneg && d <= dpeakmin && d > R_dresneg && BenMark > 0)
		return -2;
	else if ((CstateFlag == -12 || CstateFlag == -2 || CstateFlag == 5 || CstateFlag == -6 || CstateFlag == 7 || CstateFlag == -5 || CstateFlag == -4) && Tdu < 0.0 && d <= R_dcapneg && d <= dpeakmin && d > R_dresneg && BenMark < 0)
		return -2;
	else if ((CstateFlag == -2 || CstateFlag == -3 || CstateFlag == 5 || CstateFlag == -6 || CstateFlag == 7 || CstateFlag == -5 || CstateFlag == -4) && Tdu < 0.0 && d <= R_dresneg && d > delUneg && d <= dpeakmin && BenMark > 0)
		return -3;
	else if ((CstateFlag == -2 || CstateFlag == -3 || CstateFlag == 5 || CstateFlag == -6 || CstateFlag == 7 || CstateFlag == -5 || CstateFlag == -4) && Tdu < 0.0 && d <= R_dresneg && d > delUneg && d <= dpeakmin && BenMark < 0)
		return -3;

	else if ((CstateFlag == -3 || CstateFlag == -30 || CstateFlag == 31) && Tdu < 0.0 && d <= delUneg && d >= Intcpt_Xaxis_neg && d <= dpeakmin && BenMark > 0)
		return -30;
	else if ((CstateFlag == -3 || CstateFlag == -30 || CstateFlag == 31) && Tdu < 0.0 && d <= delUneg && d >= Intcpt_Xaxis_neg && d <= dpeakmin && BenMark < 0)
		return -30;

	else if ((CstateFlag == -30 || CstateFlag == -40) && Tdu < 0.0 && d < Intcpt_Xaxis_neg && d <= dpeakmin && BenMark > 0)
		return -40;
	else if ((CstateFlag == -30 || CstateFlag == -40) && Tdu < 0.0 && d < Intcpt_Xaxis_neg && d <= dpeakmin && BenMark < 0)
		return -40;
	else if ((CstateFlag == -30 || CstateFlag == 31) && (Tdu > 0.0 || Tdu < 0.0) && (d > dpeakmin || d < dpeakmax))
		return 31;
	else if ((CstateFlag == -40 || CstateFlag == 41) && (Tdu > 0.0 || Tdu < 0.0))
		return 41;

	else if ((CstateFlag == -12 || CstateFlag == -2 || CstateFlag == -3) && Tdu > 0.0)
		return -4;
	else if ((CstateFlag == -4 || CstateFlag == -5) && Tdu > 0.0 && d < dpeakmax)
		return -5;
	else if ((CstateFlag == -4 || CstateFlag == -5 || CstateFlag == -6) && Tdu < 0.0 && d > dpeakmin)
		return -6;
	else if ((CstateFlag == -6 || CstateFlag == -7 || CstateFlag == 7) && Tdu > 0.0 && d < dpeakmax)
		return -7;
	else
		return 999;

}

//Defining the peaks of the splin curve
void
GMG_CyclicReinforcedConcrete::define_peak(void)
{
	double ffmax_effective, ffmin_effective;
	double delta_R_fcappos, delta_R_fcapneg;
	double delta_R_dcappos, delta_R_dcapneg;
	//---------------Part 2: Updating the input variable for "Splineparam" function-------------------

	/* ************************************************************ **
	 * ************************************************************ **

	          Post Yielding (before capping) Damage Mode 
	               First Part, Backbone to Backbone

	 * ************************************************************ **
	 * ************************************************************ **/

	if (CstateFlag == 12 && TstateFlag == 4) {

		dpeakmax = Cstrain;
		ffmax = Cstress;

		if (flagdmg_Hardening_strength == 1) {

			delta_R_fcapneg = dmgSneg * fabs(fcapneg);
			R_dcapneg = CR_dcapneg + (fabs(Kdegneg) / (fabs(slope_neg) + fabs(Kdegneg))) * (fabs(delta_R_fcapneg) / fabs(Kdegneg));  //Updating dcapneg
			R_fcapneg = -(slope_neg * fabs(R_dcapneg) - Intcpt_slope_neg);
			R_dresneg = CR_dresneg + fabs(delta_R_fcapneg) / fabs(Kdegneg); //Updating dresneg
			Intcpt_deg_neg = CIntcpt_deg_neg - fabs(delta_R_fcapneg);
		}

		if (flagdmg_Hardening == 1) {

			dpealmax_bench = dpeakmax * (1 + delta_ratio_max_hard);

			if (dpeakmin > dpeakmin_bench) {
				dpeakmin = fmax((1 + alpha_neg) * dpeakmin, dpeakmin_bench);
				ffmin = -(slope_neg * fabs(dpeakmin) - Intcpt_slope_neg);
				if (dpeakmin < R_dcapneg) {
					ffmin = -(Kdegneg * fabs(dpeakmin) + Intcpt_deg_neg);
				}
			}
			else if (dpeakmin <= dpeakmin_bench) {
				dpeakmin = dpeakmin_bench;
				ffmin = -(slope_neg * fabs(dpeakmin) - Intcpt_slope_neg);
				if (dpeakmin < R_dcapneg) {
					ffmin = -(Kdegneg * fabs(dpeakmin) + Intcpt_deg_neg);
				}
			}
		}

	}

	if (CstateFlag == -12 && TstateFlag == -4) {
		dpeakmin = Cstrain;
		ffmin = Cstress;

		if (flagdmg_Hardening_strength == 1) {

			delta_R_fcappos = dmgSpos * fcappos;
			R_dcappos = CR_dcappos - fabs(Kdegpos) / (fabs(slope_pos) + fabs(Kdegpos))*(fabs(delta_R_fcappos) / fabs(Kdegpos));
			R_fcappos = slope_pos * fabs(R_dcappos) + Intcpt_slope_pos;
			R_drespos = CR_drespos - fabs(delta_R_fcappos) / fabs(Kdegpos);
			Intcpt_deg_pos = CIntcpt_deg_pos - fabs(delta_R_fcappos);
		}

		if (flagdmg_Hardening == 1) {
			dpeakmin_bench = dpeakmin * (1 + delta_ratio_max_hard);

			if ((dpeakmax) < dpealmax_bench) {
				dpeakmax = fmin((1 + alpha_pos) * dpeakmax, dpealmax_bench);
				ffmax = slope_pos * dpeakmax + Intcpt_slope_pos;
				if (dpeakmax > R_dcappos) {
					ffmax = Kdegpos * fabs(dpeakmax) + Intcpt_deg_pos;
				}
			}
			else if ((dpeakmax) >= dpealmax_bench) {
				dpeakmax = dpealmax_bench;
				ffmax = slope_pos * dpeakmax + Intcpt_slope_pos;
				if (dpeakmax > R_dcappos) {
					ffmax = Kdegpos * fabs(dpeakmax) + Intcpt_deg_pos;
				}
			}
		}

	}

	/* ************************************************************ **
	 * ************************************************************ **

			      Post Yielding before capping Damage Mode 
			            Second Part, inner cyclic

	 * ************************************************************ **
	 * ************************************************************ **/

	if (flagdmg_Hardening_strength == 1) {
		if ((CstateFlag == 5 && TstateFlag == 6) || (CstateFlag == -6 && TstateFlag == -7) || (CstateFlag == 7 && TstateFlag == -7)) {

			delta_R_fcappos = dmgSpos * fcappos;
			R_dcappos = CR_dcappos - fabs(Kdegpos) / (fabs(slope_pos) + fabs(Kdegpos))*(fabs(delta_R_fcappos) / fabs(Kdegpos));
			R_fcappos = slope_pos * fabs(R_dcappos) + Intcpt_slope_pos;
			R_drespos = CR_drespos - fabs(delta_R_fcappos) / fabs(Kdegpos);
			Intcpt_deg_pos = CIntcpt_deg_pos - fabs(delta_R_fcappos);
			if (dpeakmax > R_dcappos) {
				ffmax = Kdegpos * fabs(dpeakmax) + Intcpt_deg_pos;
			}
		}
	}


	if (flagdmg_Hardening_strength == 1) {
		if ((CstateFlag == -5 && TstateFlag == -6) || (CstateFlag == 6 && TstateFlag == 7) || (CstateFlag == -7 && TstateFlag == 7)) {

			delta_R_fcapneg = dmgSneg * fabs(fcapneg);
			R_dcapneg = CR_dcapneg + (fabs(Kdegneg) / (fabs(slope_neg) + fabs(Kdegneg))) * (fabs(delta_R_fcapneg) / fabs(Kdegneg));  //Updating dcapneg
			R_fcapneg = -(slope_neg * fabs(R_dcapneg) - Intcpt_slope_neg);
			R_dresneg = CR_dresneg + fabs(delta_R_fcapneg) / fabs(Kdegneg); //Updating dresneg
			Intcpt_deg_neg = CIntcpt_deg_neg - fabs(delta_R_fcapneg);
			if (dpeakmin < R_dcapneg) {
				ffmin = -(Kdegneg * fabs(dpeakmin) + Intcpt_deg_neg);
			}
		}
	}

	if (flagdmg_Hardening == 1 && (dpeakmax < dpealmax_bench)) {
		if ((CstateFlag == 5 && TstateFlag == 6) || (CstateFlag == -6 && TstateFlag == -7) || (CstateFlag == 7 && TstateFlag == -7)) {
			dpeakmax = fmin((1 + alpha_pos) * dpeakmax, dpealmax_bench);
			ffmax = slope_pos * dpeakmax + Intcpt_slope_pos;
			if (dpeakmax > R_dcappos) {
				ffmax = Kdegpos * fabs(dpeakmax) + Intcpt_deg_pos;
			}
		}
	}

	if (flagdmg_Hardening == 1 && (dpeakmin > dpeakmin_bench)) {
		if ((CstateFlag == -5 && TstateFlag == -6) || (CstateFlag == 6 && TstateFlag == 7) || (CstateFlag == -7 && TstateFlag == 7)) {
			dpeakmin = fmax((1 + alpha_neg) * dpeakmin, dpeakmin_bench);
			ffmin = -(slope_neg * fabs(dpeakmin) - Intcpt_slope_neg);
			if (dpeakmin < R_dcapneg) {
				ffmin = -(Kdegneg * fabs(dpeakmin) + Intcpt_deg_neg);
			}
		}
	}

	/* ************************************************************ **
	 * ************************************************************ **

						 End of Post Yielding Damage Mode

	 * ************************************************************ **
	 * ************************************************************ **/



	 /* ************************************************************ **
	  * ************************************************************ **

				  Post Capping (before residual) Damage Mode
				      First Part, Backbone to Backbone

	  * ************************************************************ **
	  * ************************************************************ **/

	if (CstateFlag == 2 && TstateFlag == 4) {
		dpeakmax = Cstrain;
		ffmax = Cstress;

		ratio = (dpeakmax - dcappos) / (drespos - dcappos);
		dpeakmin = (dresneg - dcapneg) * ratio + dcapneg;

		ffmin_effective = -(Kdegneg * fabs(dpeakmin) + Intcpt_deg_neg);
		ffmin = ffmin_effective;
		if (dmgSneg > 0.0) {

			delta_R_fcapneg = dmgSneg * fabs(fcapneg);
			ffmin = fmin((ffmin_effective + delta_R_fcapneg), R_fresneg);
			R_dcapneg = CR_dcapneg + (fabs(Kdegneg) / (fabs(slope_neg) + fabs(Kdegneg))) * (fabs(delta_R_fcapneg) / fabs(Kdegneg));  //Updating dcapneg
			R_dresneg = CR_dresneg + fabs(delta_R_fcapneg) / fabs(Kdegneg); //Updating dresneg
			Intcpt_deg_neg = CIntcpt_deg_neg - fabs(delta_R_fcapneg);
		}

	}
	if (CstateFlag == -2 && TstateFlag == -4) {
		dpeakmin = Cstrain;
		ffmin = Cstress;

		ratio = (dpeakmin - dcapneg) / (dresneg - dcapneg);
		dpeakmax = (drespos - dcappos) * ratio + dcappos;

		ffmax_effective = Kdegpos * fabs(dpeakmax) + Intcpt_deg_pos;
		ffmax = ffmax_effective;
		if (dmgSpos > 0.0) {

			delta_R_fcappos = dmgSpos * fcappos;
			ffmax = fmax(ffmax_effective - delta_R_fcappos, frespos);
			R_drespos = CR_drespos - fabs(delta_R_fcappos) / fabs(Kdegpos);
			R_dcappos = CR_dcappos - fabs(Kdegpos) / (fabs(slope_pos) + fabs(Kdegpos))*(fabs(delta_R_fcappos) / fabs(Kdegpos));
			Intcpt_deg_pos = CIntcpt_deg_pos - fabs(delta_R_fcappos);
		}
	}

	/* ************************************************************ **
	 * ************************************************************ **

			 Post Capping (before residual) Damage Mode
				    Second Part, Inner cyclic

	 * ************************************************************ **
	 * ************************************************************ **/

	if (flagdmg == 1) {
		if ((CstateFlag == 5 && TstateFlag == 6) || (CstateFlag == -6 && TstateFlag == -7) || (CstateFlag == 7 && TstateFlag == -7)) {

			ffmax_effective = Kdegpos * fabs(dpeakmax) + Intcpt_deg_pos;
			ffmax = ffmax_effective;
			if (dmgSpos > 0.0) {

				delta_R_fcappos = dmgSpos * fcappos;
				ffmax = fmax(ffmax_effective - delta_R_fcappos, frespos);
				R_drespos = CR_drespos - fabs(delta_R_fcappos) / fabs(Kdegpos);
				R_dcappos = CR_dcappos - fabs(Kdegpos) / (fabs(slope_pos) + fabs(Kdegpos))*(fabs(delta_R_fcappos) / fabs(Kdegpos));
				Intcpt_deg_pos = CIntcpt_deg_pos - fabs(delta_R_fcappos);
			}
		}
	}


	if (flagdmg == 1) {
		if ((CstateFlag == -5 && TstateFlag == -6) || (CstateFlag == 6 && TstateFlag == 7) || (CstateFlag == -7 && TstateFlag == 7)) {

			ffmin_effective = -(Kdegneg*fabs(dpeakmin) + Intcpt_deg_neg);
			ffmin = ffmin_effective;
			if (dmgSneg > 0.0) {

				delta_R_fcapneg = dmgSneg * fabs(fcapneg);
				ffmin = fmin((ffmin_effective + delta_R_fcapneg), R_fresneg);
				R_dcapneg = CR_dcapneg + (fabs(Kdegneg) / (fabs(slope_neg) + fabs(Kdegneg))) * (fabs(delta_R_fcapneg) / fabs(Kdegneg));  //Updating dcapneg
				R_dresneg = CR_dresneg + fabs(delta_R_fcapneg) / fabs(Kdegneg); //Updating dresneg
				Intcpt_deg_neg = CIntcpt_deg_neg - fabs(delta_R_fcapneg);
			}

		}
	}

	/* ************************************************************ **
	 * ************************************************************ **

						End of Post Capping Damage Mode

	 * ************************************************************ **
	 * ************************************************************ **/

	/* ************************************************************ **
	 * ************************************************************ **

					  Post Capping (After residual)

	 * ************************************************************ **
	 * ************************************************************ **/

	if (CstateFlag == 3 && TstateFlag == 4) {
		dpeakmax = Cstrain;
		ffmax = Cstress;
		ratio = (dpeakmax - R_drespos) / (delUpos - R_drespos);
		dpeakmin = (delUneg - R_dresneg) * ratio + R_dresneg;
		ffmin = fresneg;
	}

	if (CstateFlag == -3 && TstateFlag == -4) {
		dpeakmin = Cstrain;
		ffmin = Cstress;
		ratio = (dpeakmin - R_dresneg) / (delUneg - R_dresneg);
		dpeakmax = (delUpos - R_drespos) * ratio + R_drespos;
		ffmax = frespos;
	}


	if (CstateFlag == 30 && TstateFlag == 31) {
		dpeakmax = Cstrain;
		ffmax = Cstress;
		ratio = (dpeakmax - delUpos) / (Intcpt_Xaxis_pos - delUpos);
		dpeakmin = (Intcpt_Xaxis_neg - delUneg) * ratio + delUneg;
		ffmin = -(Kdegneg*fabs(dpeakmin) + Intcpt_res_neg);
	}

	if (CstateFlag == -30 && TstateFlag == 31) {
		dpeakmin = Cstrain;
		ffmin = Cstress;
		ratio = (dpeakmin - delUneg) / (Intcpt_Xaxis_neg - delUneg);
		dpeakmax = (Intcpt_Xaxis_pos - delUpos) * ratio + delUpos;
		ffmax = Kdegpos * fabs(dpeakmax) + Intcpt_res_pos;
	}
}
//-----------------------------My subrotin----------------------------------

// Spline parameter subroutine
void
GMG_CyclicReinforcedConcrete::splineparam(double MtoRref, double dpeakpos, double dcappose, double dpeakneg, double dcapneg)
{
	if ((fabs(dpeakpos) < fabs(dcappose)) && (fabs(dpeakneg) < fabs(dcapneg))) {
		ER = fmin(alpha_Er_Post_Yielding * pow((ffmax - ffmin) / (dpeakmax - dpeakmin) / ((Kepos + Keneg) / 2), beta_Er_Post_Yielding), Er_Post_Yielding);
		KrR = Kr_Post_Yielding;
		KuR = Kun_Post_Yielding;
	}
	else {
		ER = fmin(alpha_Er_Post_Capping * pow((ffmax - ffmin) / (dpeakmax - dpeakmin) / ((Kepos + Keneg) / 2), beta_Er_Post_Capping), Er_Post_Capping);
		KrR = Kr_Post_Capping;
		KuR = Kun_Post_Capping;
	}
}

void
GMG_CyclicReinforcedConcrete::spline_curve(double dirtag, double P0x, double P0y, double P4x, double P4y, double Krel, double Kun, double E, double x)
{
	double delx, C, P1x, P2x, P3x, P1y, P2y, P3y;
	double a, b, c, Q, R, theta, u;
	double B11, B21, B31, B41, B51;
	double B12, B22, B32, B42, B52;
	const double PI = 3.141592653589793238463;

	if (dirtag == 0) { // NP dir
		x = dpeakmax - (x - dpeakmin);
	}

	delx = (P4x - P0x) / 4;
	if (P4x * 100 > 4) {
		C = 0.1;
	}
	else {
		C = 0.1;
	}
	P1x = P0x + C / 4 * (P4x - P0x);
	P2x = P0x + 2 * delx;
	P3x = P4x - C / 4 * (P4x - P0x);
	P3y = P4y - Kun * (C*delx);
	P1y = P0y + Krel * (C*delx);
	P2y = 0.5 * C*delx*(2.0 * (Kun - Krel)) - 0.5 * E / delx + .5 * (P0y + P4y);

	// mapping x to u
	a = (P0x - P4x)*(9.0 / 2.0 * C - 3.0) / ((P0x - P4x)*(2.0 - 3.0 * C));
	b = 3.0 / 2.0 * C*(P4x - P0x) / ((P0x - P4x)*(2.0 - 3.0 * C));
	c = (P0x - x) / ((P0x - P4x)*(2.0 - 3.0 * C));
	Q = (pow(a, 2.0) - 3.0 * b) / 9.0;
	R = (2.0 * pow(a, 3.0) - 9.0 * a*b + 27.0 * c) / 54.0;
	theta = acos(R / sqrt(pow(Q, 3.0)));
	u = -2.0 * sqrt(Q)*cos((theta - 2.0 * PI) / 3.0) - a / 3.0;
	// solve basis function and its derivative
	if (u < 0.5) {
		B11 = -8.0 * pow((u - 0.5), 3.0);
		B21 = 14.0 * u*(0.428571 - 1.28571*u + pow(u, 2.0));
		B31 = -8.0 * pow(u, 2.0) * (u - 0.75);
		B41 = 2.0 * pow(u, 3.0);
		B51 = 0.0;
		B12 = -24.0 * pow((u - 0.5), 2.0);
		B22 = 6.0 * (1.0 - 6.0 * u + 7.0 * pow(u, 2.0));
		B32 = 12.0 * u*(1.0 - 2.0 * u);
		B42 = 6.0 * pow(u, 2.0);
		B52 = 0.0;
	}
	else {
		B11 = 0.0;
		B21 = -2.0 * pow((u - 1.0), 3.0);
		B31 = 8.0 * pow((1.0 - u), 2.0) * (u - 0.25);
		B41 = 2.0 - 12.0 * u + 24.0 * pow(u, 2.0) - 14.0 * pow(u, 3.0);
		B51 = 8.0 * pow((u - 0.5), 3.0);
		B12 = 0.0;
		B22 = -6.0 * pow((1.0 - u), 2.0);
		B32 = 12.0 * (1.0 - 3.0 * u + 2.0 * pow(u, 2.0));
		B42 = -12.0 + 48.0 * u - 42.0 * pow(u, 2.0);
		B52 = 24.0 * pow((u - 0.5), 2.0);
	}
	fspl = P0y * B11 + P1y * B21 + P2y * B31 + P3y * B41 + P4y * B51;
	fderspl = P0y * B12 + P1y * B22 + P2y * B32 + P3y * B42 + P4y * B52;
	dderspl = P0x * B12 + P1x * B22 + P2x * B32 + P3x * B42 + P4x * B52;
	ek = fderspl / dderspl;
	if (dirtag == 0) { // NP dir
		fspl = -fspl + P0y + P4y;
	}
}

void
GMG_CyclicReinforcedConcrete::checkEnvelope(void)
{

	if (f >= 0.0 && d >= 0.0) {

		if (d >= dpeakmax && d <= R_dcappos)
		{
			TstateFlag = 12;
			ek = slope_pos;
			d = dpeakmax;
			f = ffmax;
		}

		if (d >= dpeakmax && d >= R_dcappos && d <= R_drespos)
		{
			TstateFlag = 2;
			ek = Kdegpos;
			d = dpeakmax;
			f = ffmax;
		}

		if (d >= R_dcappos && d <= R_drespos && d < dpeakmax && f >= Kdegpos * fabs(d) + Intcpt_deg_pos) {
			TstateFlag = 2;
			ek = Kdegpos;
			f = Kdegpos * fabs(d) + Intcpt_deg_pos;
			dpeakmax = d;
		}

		if (d >= dpeakmax && d >= R_drespos && d < delUpos)
		{
			TstateFlag = 3;
			ek = 0.0001;
			d = dpeakmax;
			f = frespos;
		}

		if (d >= dpeakmax && d >= delUpos && d <= Intcpt_Xaxis_pos)
		{
			TstateFlag = 30;
			ek = Kdegpos;
			d = dpeakmax;
			f = ffmax;
		}
	}

	else if (f < 0.0 && d < 0.0)
	{
		if (d <= dpeakmin && d >= R_dcapneg)
		{
			TstateFlag = -12;
			ek = slope_neg;
			d = dpeakmin;
			f = ffmin;
		}
		if (d <= dpeakmin && d <= R_dcapneg && d >= R_dresneg)
		{
			TstateFlag = -2;
			ek = Kdegneg;
			d = dpeakmin;
			f = ffmin;
		}

		if (d > dpeakmin && d <= R_dcapneg && d >= R_dresneg && f <= -(Kdegneg*fabs(d) + Intcpt_deg_neg))
		{
			TstateFlag = -2;
			ek = Kdegneg;
			f = -(Kdegneg*fabs(d) + Intcpt_deg_neg);
			dpeakmin = d;
		}

		if (d <= dpeakmin && d <= R_dresneg && d > delUneg)
		{
			TstateFlag = -3;
			ek = 0.0001;
			d = dpeakmin;
			f = fresneg;
		}

		if (d <= dpeakmin && d <= delUneg && d >= Intcpt_Xaxis_neg)
		{
			TstateFlag = -30;
			ek = Kdegneg;
			d = dpeakmin;
			f = ffmin;
		}
	}
}
//************************************************************************************************************

//**********************************************My subrotin***************************************************

//This section has been added to update the damage after entering case 2 or -2 for the first time
void
GMG_CyclicReinforcedConcrete::update_damage(void)
{
	double Energy_Coe;
	spos = (f + fabs(f)) / 2;
	sneg = (f - fabs(f)) / 2;
	if (d >= dpeakmax || d <= dpeakmin) {
		Energy_Coe = C3;
	}
	else if (Tdu <= 0 && f >= 0) {
		Energy_Coe = C1;
	}
	else if (Tdu <= 0 && f < 0) {
		Energy_Coe = C2;
	}
	else if (Tdu >= 0 && f >= 0) {
		Energy_Coe = C2;
	}
	else if (Tdu >= 0 && f < 0) {
		Energy_Coe = C1;
	}
	Et = CEt + Energy_Coe * fabs(Cspos + spos) * fabs(Tdu) / 2;
	Ec = CEc + Energy_Coe * fabs(Csneg + sneg) * fabs(Tdu) / 2;
	Epos = Et + Ec;
	Eneg = Et + Ec;

	if (Tdu < 0) {

		if ((CstateFlag == -5 && TstateFlag == -6) || (CstateFlag == 6 && TstateFlag == 7) || (CstateFlag == -7 && TstateFlag == 7) || (CstateFlag == 2 && TstateFlag == 4) || (CstateFlag == 12 && TstateFlag == 4)) {

			spos_neg = Cspos_neg = sneg_neg = Csneg_neg = Et_neg = CEt_neg = Ec_neg = CEc_neg = Enegnorm = dmgSneg = 0.0;

		}

		if (Epos >= Ed0pos) {

			spos_pos = (f + fabs(f)) / 2;
			sneg_pos = (f - fabs(f)) / 2;
			if (d >= dpeakmax || d <= dpeakmin) {
				Energy_Coe = C3;
			}
			else if (Tdu <= 0 && f >= 0) {
				Energy_Coe = C1;
			}
			else if (Tdu <= 0 && f < 0) {
				Energy_Coe = C2;
			}
			Et_pos = CEt_pos + Energy_Coe * fabs(Cspos_pos + spos_pos) * fabs(Tdu) / 2;
			Ec_pos = CEc_pos + Energy_Coe * fabs(Csneg_pos + sneg_pos) * fabs(Tdu) / 2;
			Epos_pos = Et_pos + Ec_pos;
			Eposnorm = fmax(0.0, (Epos_pos) / (Ed1pos - Ed0pos));


			if (flagdmg == 1) {
				dmgSpos = fmin(solpe_post_capping*Eposnorm, 1); //cdf(d, Eposnorm);
			}
			if (flagdmg_Hardening_strength == 1) {
				dmgSpos = fmin(solpe_post_yielding*Eposnorm, 1); //cdf(d, Eposnorm);
			}
		}
	}

	if (Tdu > 0) {

		if ((CstateFlag == 5 && TstateFlag == 6) || (CstateFlag == -6 && TstateFlag == -7) || (CstateFlag == 7 && TstateFlag == -7) || (CstateFlag == -2 && TstateFlag == -4) || (CstateFlag == -12 && TstateFlag == -4)) {

			spos_pos = Cspos_pos = sneg_pos = Csneg_pos = Et_pos = CEt_pos = Ec_pos = CEc_pos = Eposnorm = dmgSpos = 0.0;

		}

		if (Epos >= Ed0pos) {

			spos_neg = (f + fabs(f)) / 2;
			sneg_neg = (f - fabs(f)) / 2;
			if (d >= dpeakmax || d <= dpeakmin) {
				Energy_Coe = C3;
			}
			else if (Tdu >= 0 && f < 0) {
				Energy_Coe = C1;
			}
			else if (Tdu >= 0 && f >= 0) {
				Energy_Coe = C2;
			}
			Et_neg = CEt_neg + Energy_Coe * fabs(Cspos_neg + spos_neg) * fabs(Tdu) / 2;
			Ec_neg = CEc_neg + Energy_Coe * fabs(Csneg_neg + sneg_neg) * fabs(Tdu) / 2;
			Eneg_neg = Et_neg + Ec_neg;
			Enegnorm = fmax(0.0, (Eneg_neg) / (Ed1neg - Ed0neg));

			if (flagdmg == 1) {
				dmgSneg = fmin(solpe_post_capping * Enegnorm, 1); //cdf(d, Enegnorm);
			}
			if (flagdmg_Hardening_strength == 1) {
				dmgSneg = fmin(solpe_post_yielding * Enegnorm, 1); //cdf(d, Enegnorm);
			}
		}
	}
}
//************************************************************************************************************


void
GMG_CyclicReinforcedConcrete::update_damage_hardeingin(void)
{
	double Energy_Coe;
	spos_hard = (f + fabs(f)) / 2;
	sneg_hard = (f - fabs(f)) / 2;
	Energy_Coe = 1.0;
	Et_hard = CEt_hard + Energy_Coe * fabs(Cspos_hard + spos_hard) * fabs(Tdu) / 2;
	Ec_hard = CEc_hard + Energy_Coe * fabs(Csneg_hard + sneg_hard) * fabs(Tdu) / 2;
	Epos_hard = Et_hard + Ec_hard;
	Eneg_hard = Et_hard + Ec_hard;

	if (Tdu < 0 && dpeakmax > dypos) {

		if ((CstateFlag == -5 && TstateFlag == -6) || (CstateFlag == 6 && TstateFlag == 7) || (CstateFlag == -7 && TstateFlag == 7) || (CstateFlag == 12 && TstateFlag == 4)) {

			spos_neg_hard = Cspos_neg_hard = sneg_neg_hard = Csneg_neg_hard = Et_neg_hard = CEt_neg_hard = Ec_neg_hard = CEc_neg_hard = 0.0;
			Et_neg_hard = CEt_neg_hard = Ec_neg_hard = CEc_neg_hard = Eneg_neg_hard = Enegnorm_hard = 0.0;
		}

		if (Epos_hard >= Ed0pos_hard) {

			spos_pos_hard = (f + fabs(f)) / 2;
			sneg_pos_hard = (f - fabs(f)) / 2;
			Energy_Coe = 1.0;

			Et_pos_hard = CEt_pos_hard + Energy_Coe * fabs(Cspos_pos_hard + spos_pos_hard) * fabs(Tdu) / 2;
			Ec_pos_hard = CEc_pos_hard + Energy_Coe * fabs(Csneg_pos_hard + sneg_pos_hard) * fabs(Tdu) / 2;
			Epos_pos_hard = Et_pos_hard + Ec_pos_hard;
			Eposnorm_hard = fmax(0.0, (Epos_pos_hard) / (Ed1pos_hard - Ed0pos_hard));
			delta_pos_hard = solpe_post_capping * Eposnorm_hard /* + Cdelta_pos_hard*/;
			alpha_pos = delta_ratio_max_hard * Eposnorm_hard;
		}
	}

	if (Tdu > 0 && dpeakmin < dyneg) {

		if ((CstateFlag == 5 && TstateFlag == 6) || (CstateFlag == -6 && TstateFlag == -7) || (CstateFlag == 7 && TstateFlag == -7) || (CstateFlag == -12 && TstateFlag == -4)) {

			spos_pos_hard = Cspos_pos_hard = sneg_pos_hard = Csneg_pos_hard = Et_pos_hard = CEt_pos_hard = Ec_pos_hard = CEc_pos_hard = 0.0;
			Et_pos_hard = CEt_pos_hard = Ec_pos_hard = CEc_pos_hard = Epos_pos_hard = Eposnorm_hard = 0.0;
		}

		if (Epos_hard >= Ed0pos_hard) {

			spos_neg_hard = (f + fabs(f)) / 2;
			sneg_neg_hard = (f - fabs(f)) / 2;
			Energy_Coe = 1.0;

			Et_neg_hard = CEt_neg_hard + Energy_Coe * fabs(Cspos_neg_hard + spos_neg_hard) * fabs(Tdu) / 2;
			Ec_neg_hard = CEc_neg_hard + Energy_Coe * fabs(Csneg_neg_hard + sneg_neg_hard) * fabs(Tdu) / 2;
			Eneg_neg_hard = Et_neg_hard + Ec_neg_hard;
			Enegnorm_hard = fmax(0.0, (Eneg_neg_hard) / (Ed1neg_hard - Ed0neg_hard));
			delta_neg_hard = solpe_post_capping * Enegnorm_hard /* + Cdelta_neg_hard */;

			alpha_neg = delta_ratio_max_hard * Enegnorm_hard;

		}
	}
}
