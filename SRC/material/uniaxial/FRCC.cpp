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

//This code is developed and written by FERAS L. KHLEF (feraskhlef@gmail.com)
// $Version: 1.0 $
// $Date: 2020-03-03$
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FRCC.h $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FRCC.cpp $

// Description: This file contains the class implementation for FRCC. 

#include <FRCC.h>
#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <iostream>
#include <fstream>
#include <elementAPI.h>
using namespace std;


void* OPS_FRCC()
{
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 25) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: uniaxialMaterial FRCC tag? EC? ET1? FT1? ET2? FT2? ETU? FTU? STU? PARAT1? PARAT2? PARAT3? PARAT4? PARAT5? PARAT6?";
		opserr << "ECP? FCP? ECU? RC? PARAC1? PARAC2? PARAC3? PARAC4? PARAC5? PARAC6?\n";
		return 0;
	}

	int tag;
	numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
		return 0;
	}

	double data[24];
	numdata = 24;
	if (OPS_GetDoubleInput(&numdata, data)) {
		return 0;
	}

	UniaxialMaterial* mat = new FRCC(tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15], data[16], data[17], data[18], data[19], data[20], data[21], data[22], data[23]);
	if (mat == 0) {
		opserr << "WARNING: failed to create FRCC material\n";
		return 0;
	}

	return mat;
}
FRCC::FRCC
(int tag, double EC, double ET1, double FT1, double ET2, double FT2, double ETU, double FTU, double STU, double PARAT1, double PARAT2, double PARAT3,
	double PARAT4, double PARAT5, double PARAT6, double ECP, double FCP, double ECU, double PARAC1, double PARAC2, double PARAC3, double PARAC4,
	double PARAC5, double PARAC6, double RC)
	:UniaxialMaterial(tag, MAT_TAG_FRCC),
	Ec(EC), et1(ET1), ft1(FT1), et2(ET2), ft2(FT2), etu(ETU), ftu(FTU), stu(STU), paraT1(PARAT1), paraT2(PARAT2), paraT3(PARAT3), paraT4(PARAT4),
	paraT5(PARAT5), paraT6(PARAT6), ecp(ECP), fcp(FCP), ecu(ECU), paraC1(PARAC1), paraC2(PARAC2), paraC3(PARAC3), paraC4(PARAC4),
	paraC5(PARAC5), paraC6(PARAC6), rc(RC)
{
	//------------------------------------------------------
	// Set all history and state variables to initial values
	//------------------------------------------------------
	this->revertToStart();

	//--------------------------------
	//calculations for the compression
	//--------------------------------
	nc = Ec * ecp / fcp;

	Dfunc(rc, ecu / ecp, nc);

	Ycr = nc * (ecu / ecp) / D;
	Zcr = (1.0 - pow(ecu / ecp, rc)) / (pow(D, 2.0));
	Xsp = ecu / ecp - Ycr / nc / Zcr;

	// These steps are necessary to ensure that the following compression inputs are always inputed negative
	if (fcp > 0.0) {
		fcp = -fcp;
	}
	if (ecp > 0.0) {
		ecp = -ecp;
	}
	if (ecu > 0.0) {
		ecu = -ecu;
	}
	//----------------------------
	//calculations for the tension
	//----------------------------
	Eh = (ft2 - ft1) / (et2 - et1);  // the slope of the hardening branch
	EsecDesce = (ftu - ft2) / (etu - et2);
	Xcrk = etu + (ftu / stu);
	RRR = (-stu - EsecDesce) / (EsecDesce - Ec);
	if (RRR < 0.0) {
		RRR = 0.00001;
	}
	if ( RRR > 300) {
	RRR = 300.0;
	}
	AAA = (EsecDesce - Ec) / pow(fabs(etu - et2), RRR);

}
FRCC::FRCC() :UniaxialMaterial(0, MAT_TAG_FRCC),
Ec(0.0), et1(0.0), ft1(0.0), et2(0.0), ft2(0.0), etu(0.0), ftu(0.0), stu(0.0), paraT1(0.0), paraT2(0.0), paraT3(0.0), paraT4(0.0),
paraT5(0.0), paraT6(0.0), ecp(0.0), fcp(0.0), ecu(0.0), paraC1(0.0), paraC2(0.0), paraC3(0.0), paraC4(0.0), paraC5(0.0), paraC6(0.0), rc(0.0),
Cepstul(0.0), CsigtMax(0.0), Cepscul(0.0), CsigcMax(0.0), CepstMax(0.0), CepscMax(0.0), CepstMaxStar(0.0), CepscMaxStar(0.0), CsigtMaxStar(0.0),
CsigcMaxStar(0.0), CepstulStar(0.0), CepsculStar(0.0), CsigtulStar(0.0), CsigculStar(0.0), CFlagT(0), CFlagC(0), CepsculStar1(0.0),
CsigculStar1(0.0), CepsculStar2(0.0), CsigculStar2(0.0), CepsculStar3(0.0), CsigculStar3(0.0), CepstulStar1(0.0), CsigtulStar1(0.0),
CepstulStar2(0.0), CsigtulStar2(0.0), CepstulStar3(0.0), CsigtulStar3(0.0), CTI(0.0), CepstRe(0.0), CepscRe(0.0), Cepst0(0.0),
Cepsc0(0.0), Cid1(0), Cid2(0), Cid3(0), Cid4(0), Cid2_1(0), Cid4_1(0)
{
	// Set trial values
	this->revertToLastCommit();
}

FRCC::~FRCC()
{
	// Does nothing// No dynamic variables are used so a destructor is not required.
}


int FRCC::setTrialStrain(double strain, double strainRate)
{
	// Set trial strain
	Tepstul = Cepstul;
	TsigtMax = CsigtMax;
	Tepscul = Cepscul;
	TsigcMax = CsigcMax;
	TepstMax = CepstMax;
	TepscMax = CepscMax;
	TepstMaxStar = CepstMaxStar;
	TepscMaxStar = CepscMaxStar;
	TsigtMaxStar = CsigtMaxStar;
	TsigcMaxStar = CsigcMaxStar;
	TepstulStar = CepstulStar;
	TepsculStar = CepsculStar;
	TsigtulStar = CsigtulStar;
	TsigculStar = CsigculStar;
	TepsculStar1 = CepsculStar1;
	TsigculStar1 = CsigculStar1;
	TepsculStar2 = CepsculStar2;
	TsigculStar2 = CsigculStar2;
	TepsculStar3 = CepsculStar3;
	TsigculStar3 = CsigculStar3;
	TepstulStar1 = CepstulStar1;
	TsigtulStar1 = CsigtulStar1;
	TepstulStar2 = CepstulStar2;
	TsigtulStar2 = CsigtulStar2;
	TepstulStar3 = CepstulStar3;
	TsigtulStar3 = CsigtulStar3;
	TepstRe = CepstRe;
	TepscRe = CepscRe;
	Tepst0 = Cepst0;
	Tepsc0 = Cepsc0;
	TTI = CTI;
	TFlagT = CFlagT;
	TFlagC = CFlagC;
	Tid1 = Cid1;
	Tid2 = Cid2;
	Tid3 = Cid3;
	Tid4 = Cid4;
	Tid2_1 = Cid2_1;
	Tid4_1 = Cid4_1;
	 
	//Set trial strain
	//-----------------
	Tstrain = strain;
	
	double DeltaStrain = Tstrain - Cstrain;

	// Skip trivial strain increments
	//-----------------
	if (fabs(DeltaStrain) < DBL_EPSILON)
	{
		Tstress = Cstress;
		Ttangent = Ctangent;
		return 0;
	}

	//-------------------
	// Code starts here
	//-------------------
	if (Tstrain > Xcrk) { // FRCC passed tension limits
		FlagT0 = 1; // FRCC failed in tension
		FlagC0 = 0; // FRCC not failed yet in compression
	}
	if (Tstrain < (Xsp * ecp)) { // FRCC passed compression limits
		FlagT0 = 1; // FRCC failed in tension
		FlagC0 = 1; // FRCC failed in compression
	}
	if (FlagC0 == 1) { // if failed in compression then no stress capacity
		Tstress = 0.0;
		Ttangent = 0.0;
	}
	else {
		if (Tstrain < 0.0) { // Under compression
			if (Tstrain > Cstrain) { //unloading in compression 
				Tid4 = 0;
				Tid4_1 = 0;
				TFlagC = 1;

				RegExpr3(TepscMax, ecu, Ec, paraC3);
				TepsculStar = TepscMaxStar - TsigcMaxStar / Esec;
				if (TepsculStar >= 0.0) {
					TepsculStar = 0.0;
					Tid1 = 0;
					Tid2 = 1;
				}
				else {
					Tid2 = 0;
				}
				if (Tstrain <= TepsculStar) { // unloading from (TepscMax, TsigcMax) toward (epstul,0)
					RegExpr4(TepscMax, ecu, Ec, paraC4);
					if (FlagT0 == 1) {
						Eul = 0.0;
					}
					TstressAndTtangent(TepscMaxStar, TsigcMaxStar, Ec, TepsculStar, 0.0, Eul, Tstrain);

					TepsculStar1 = Tstrain;
					TsigculStar1 = Tstress;
				}
				else {// This is a gradual crack closure. unloading from (epscul,0) to either a point on the tension envelope or to (TepstMax, TsigtMax)
					if (FlagT0 == 1) {
						Tid2 = 1;
						Tid4 = 0;
						Tstress = 0.0;
						Ttangent = 0.0;
					}
					else {
						Tid2 = 1;
						Tid4 = 0;
						if (Tid1 == 0) { //unloading from (epscul,0) to a point on the tension envelope (0,sig)

							Tepst0 = paraC5* fabs(TepscMax); // shifting strain (epst0,0) for tension side.
							RegExpr4(TepscMax, ecu, Ec, paraC4);
							EnvTen(0.0, Tepst0);

							if (Tid2_1 == 0) {
								TstressAndTtangent(TepsculStar, 0.0, Eul, 0.0, Tstress, Ttangent, Tstrain);
							}
							else {
								TstressAndTtangent(TepsculStar3, TsigculStar3, Eul, 0.0, Tstress, Ttangent, Tstrain);
							}
						}
						else { // unloading from (epscul,0) to (TepstMax, TsigtMax)

							RegExpr4(TepscMax, ecu, Ec, paraC4);

							fnewEnew(TepstMax, Tepstul, TsigtMax, paraT2);

							if (Tid2_1 == 0) {
								TstressAndTtangent(TepsculStar, 0.0, Eul, TepstMax, fnew, Enew, Tstrain);
							}
							else {
								TstressAndTtangent(TepsculStar3, TsigculStar3, Eul, TepstMax, fnew, Enew, Tstrain);
							}
						}
					}
					TepsculStar2 = Tstrain;
					TsigculStar2 = Tstress;
				}
			} //Unloading in compression
			else { // Reloading or monotonic within the compression side
				if (Tstrain >= TepscRe && TFlagC == 1) { // Reloading

					if (Tstrain >= TepscMax) {// Reloading before (TepscMax, fnew)

						if (Tstrain <= TepsculStar1 && Tid2 == 0 && Tid4 == 0) {

							fnewEnew(TepscMax, Tepscul, TsigcMax, paraC2); //  Tepstul to TepstulStar1, and Tsigtul to TsigtulStar
							TstressAndTtangent(TepsculStar1, TsigculStar1, Ec, TepscMax, fnew, Enew, Tstrain);

							TepscMaxStar = Tstrain;
							TsigcMaxStar = Tstress;
							TTI = Ttangent;
						}
						if (Tid2 == 1) {
							if (FlagT0 == 1) {
								Tid2_1 = 0;
								fnewEnew(TepscMax, Tepscul, TsigcMax, paraC2);

								TstressAndTtangent(TepsculStar2, 0.0, 0.0, TepscMax, fnew, Enew, Tstrain);

								TepscMaxStar = Tstrain;
								TsigcMaxStar = Tstress;
								TTI = Ttangent;
							}
							else {
								epsca = TepsculStar2 - TsigculStar2 / Ec;
								RegExpr6(TepscMax, ecu, Ec, paraC6);
								if (Tstrain >= epsca) {
									Tid2_1 = 1;
									TstressAndTtangent(TepsculStar2, TsigculStar2, Ec, epsca, 0.0, Ea2, Tstrain);
									TepsculStar3 = Tstrain;
									TsigculStar3 = Tstress;
								}
								else {
									Tid2_1 = 0;
									fnewEnew(TepscMax, Tepscul, TsigcMax, paraC2);
									TstressAndTtangent(epsca, 0.0, Ea2, TepscMax, fnew, Enew, Tstrain);

									TepscMaxStar = Tstrain;
									TsigcMaxStar = Tstress;
									TTI = Ttangent;
								}
							}
						}
						if (Tid4 == 1) { // These code lines complete the unloading curve from epstul to comp envelope
							RegExpr4(TepstMax, etu, Ec, paraT4);
							fnewEnew(TepscMax, Tepscul, TsigcMax, paraC2);

							if (Tid4_1 == 0) {
								if (FlagT0 == 1) {
									TstressAndTtangent(TepstMaxStar, 0.0, 0.0, TepscMax, fnew, Enew, Tstrain);
								}
								else {
									TstressAndTtangent(TepstulStar, 0.0, Eul, TepscMax, fnew, Enew, Tstrain);
								}
							}
							else {
								TstressAndTtangent(TepstulStar3, TsigtulStar3, Eul, TepscMax, fnew, Enew, Tstrain);
							}

							TepscMaxStar = Tstrain;
							TsigcMaxStar = Tstress;
							TTI = Ttangent;
						}
					}
					else { //  reloading curves after TepscMax point (from (TepscMax,fnew) to (TepscRe, TsigcRe))
						Tid2 = 0;
						Tid3 = 1;
						Tid4 = 0;
						Tid4_1 = 0;

						EnvComp(TepscRe, Tepsc0);
						fnewEnew(TepscMax, Tepscul, TsigcMax, paraC2);
						TstressAndTtangent(TepscMax, fnew, TTI, TepscRe, Tstress, Ttangent, Tstrain);

						TepscMaxStar = Tstrain;
						TsigcMaxStar = Tstress;
					}
				}
				else { // monotonic compression 
					Tid3 = 1;
					TFlagC = 0;
					EnvComp(Tstrain, Tepsc0);

					TepscMax = Tstrain;
					TsigcMax = Tstress;

					TepscMaxStar = Tstrain;
					TsigcMaxStar = Tstress;

					TepscRe = paraC1 * TepscMax;
					RegExpr3(TepscMax, ecu, Ec, paraC3);

					Tepscul = TepscMax - TsigcMax / Esec;
					if (Tepscul > 0.0) {
						Tepscul = 0.0;
					}
				}
			}
		} // Under Compression
		//------------------------------------------------------------------------------------------------------
		//------------------------------------------------------------------------------------------------------
		//------------------------------------------------------------------------------------------------------
		//------------------------------------------------------------------------------------------------------
		else { // Under tension
			if (FlagT0 == 1) { // section is failed in tension
				if (Tstrain < Cstrain) { //unloading within tension side
					Tid4 = 1;
					if (Tid3 == 0) {
						Tepsc0 = paraC5* TepstMax;
						EnvComp(0.0, Tepsc0);

						if (Tid4_1 == 0) {
							TstressAndTtangent(TepstMaxStar, 0.0, 0.0, 0.0, Tstress, Ttangent, Tstrain);
						}
						else {
							TstressAndTtangent(TepstulStar3, TsigtulStar3, 0.0, 0.0, Tstress, Ttangent, Tstrain);
						}
					}
					else {
						fnewEnew(TepscMax, Tepscul, TsigcMax, paraC2);

						if (Tid4_1 == 0) {
							TstressAndTtangent(TepstMaxStar, 0.0, 0.0, TepscMax, fnew, Enew, Tstrain);
						}
						else {
							TstressAndTtangent(TepstulStar3, TsigtulStar3, 0.0, TepscMax, fnew, Enew, Tstrain);
						}
					}
					TepstulStar2 = Tstrain;
					TsigtulStar2 = Tstress;
				}
				else {// Reloading or monotonic within tension side

					if (Tid4 == 1) {
						epsta = TepstulStar2 - TsigtulStar2 / Ec;

						if (Tstrain <= epsta) {
							Tid4_1 = 1;
							TstressAndTtangent(TepstulStar2, TsigtulStar2, Ec, epsta, 0.0, 0.0, Tstrain);
							TepstulStar3 = Tstrain;
							TsigtulStar3 = Tstress;
						}
						else {
							Tid4_1 = 0;
							Tstress = 0.0;
							Ttangent = 0.0;
							TepstMaxStar = Tstrain;
						}
					}
					else {
						Tstress = 0.0;
						Ttangent = 0.0;
						TepstMaxStar = Tstrain;
					}
				}
			}
			//---------------------------------------------------------------------------------------------------------------------------------------
			else { // section is not failed yet in tension
				if (Tstrain < Cstrain) { //unloading within tension side
					Tid2 = 0;
					Tid2_1 = 0;
					TFlagT = 1;

					RegExpr3(TepstMax, etu, Ec, paraT3);
					TepstulStar = TepstMaxStar - TsigtMaxStar / Esec;

					if (TepstulStar < 0.0) {
						TepstulStar = 0.0;
						Tid4 = 1;
					}
					else {
						Tid4 = 0;
					}
					if (Tstrain >= TepstulStar) { // unloading curve from the unloading point to epstul
						RegExpr4(TepstMax, etu, Ec, paraT4);
						TstressAndTtangent(TepstMaxStar, TsigtMaxStar, Ec, TepstulStar, 0.0, Eul, Tstrain);

						TepstulStar1 = Tstrain;
						TsigtulStar1 = Tstress;
					}
					else { // unloading from epstul to a target point on compression side
						Tid4 = 1;
						Tid2 = 0;
						
						if (Tid3 == 0) {
							Tepsc0 = paraT5 * TepstMax;
							RegExpr4(TepstMax, etu, Ec, paraT4);
							EnvComp(0.0, Tepsc0);
						
							if (Tid4_1 == 0) {
								TstressAndTtangent(TepstulStar, 0.0, Eul, 0.0, Tstress, Ttangent, Tstrain);
							}
							else {
								TstressAndTtangent(TepstulStar3, TsigtulStar3, Eul, 0.0, Tstress, Ttangent, Tstrain);
							}
						}
						else {
							RegExpr4(TepstMax, etu, Ec, paraT4);
							fnewEnew(TepscMax, Tepscul, TsigcMax, paraC2);
						
							if (Tid4_1 == 0) {
								TstressAndTtangent(TepstulStar, 0.0, Eul, TepscMax, fnew, Enew, Tstrain);
							}
							else {
								TstressAndTtangent(TepstulStar3, TsigtulStar3, Eul, TepscMax, fnew, Enew, Tstrain);
								
							}
						}
						TepstulStar2 = Tstrain;
						TsigtulStar2 = Tstress;
					}
				}
				else {// Reloading or monotonic within tension side
					
					if (Tstrain <= TepstRe && TFlagT == 1) { //Reloading 

						if (Tstrain <= TepstMax) {  // reloading curves before TepstMax point

							if (Tstrain >= TepstulStar1 && Tid2 == 0 && Tid4 == 0) {

								fnewEnew(TepstMax, Tepstul, TsigtMax, paraT2);
								TstressAndTtangent(TepstulStar1, TsigtulStar1, Ec, TepstMax, fnew, Enew, Tstrain);

								TepstMaxStar = Tstrain;
								TsigtMaxStar = Tstress;
								TTI = Ttangent;
							}
							if (Tid4 == 1) {
								epsta = TepstulStar2 - TsigtulStar2 / Ec;
								RegExpr6(TepstMax, etu, Ec, paraT6);
								if (Tstrain <= epsta) {
									Tid4_1 = 1;
									TstressAndTtangent(TepstulStar2, TsigtulStar2, Ec, epsta, 0.0, Ea2, Tstrain);
									TepstulStar3 = Tstrain;
									TsigtulStar3 = Tstress;
								}
								else {
									Tid4_1 = 0;
									fnewEnew(TepstMax, Tepstul, TsigtMax, paraT2);
									TstressAndTtangent(epsta, 0.0, Ea2, TepstMax, fnew, Enew, Tstrain);

									TepstMaxStar = Tstrain;
									TsigtMaxStar = Tstress;
									TTI = Ttangent;
								}
							}
							if (Tid2 == 1) {
								RegExpr4(TepscMax, ecu, Ec, paraC4);
								fnewEnew(TepstMax, Tepstul, TsigtMax, paraT2);

								if (Tid2_1 == 0) {
									TstressAndTtangent(TepsculStar, 0.0, Eul, TepstMax, fnew, Enew, Tstrain);
								}
								else {
									TstressAndTtangent(TepsculStar3, TsigculStar3, Eul, TepstMax, fnew, Enew, Tstrain);
								}

								TepstMaxStar = Tstrain;
								TsigtMaxStar = Tstress;
								TTI = Ttangent;
							}
						}
						else { // reloading curves after TepstMax point (from (TepstMax,fnew) to (TepstRe, sigtRe))
							Tid1 = 1;
							Tid2 = 0;
							Tid4 = 0;
							Tid2_1 = 0;
							EnvTen(TepstRe, Tepst0);
							fnewEnew(TepstMax, Tepstul, TsigtMax, paraT2);
							TstressAndTtangent(TepstMax, fnew, TTI, TepstRe, Tstress, Ttangent, Tstrain);

							TepstMaxStar = Tstrain;
							TsigtMaxStar = Tstress;
						}
					}
					else { //monotonic tension
						Tid1 = 1;
						TFlagT = 0;
						EnvTen(Tstrain, Tepst0);

						TepstMax = Tstrain;
						TsigtMax = Tstress;

						TepstMaxStar = Tstrain;
						TsigtMaxStar = Tstress;

						TepstRe = paraT1 * TepstMax;
						RegExpr3(TepstMax, etu, Ec, paraT3);

						Tepstul = TepstMax - TsigtMax / Esec;
						if (Tepstul < 0.0) {
							Tepstul = 0.0;
						}
					}
				}
			}
		}
	}
	return 0;
}
// FRCC Functions
void FRCC::RegExpr3(double epsun, double epsu, double EE, double PARA3) {
	Esec = EE * PARA3 / (fabs(epsun / epsu) + PARA3);
}
void FRCC::RegExpr4(double epsun, double epsu, double EE, double PARA4) {
	Eul = EE * PARA4 / (fabs(epsun / epsu) + PARA4);
}
void FRCC::RegExpr6(double epsun, double epsu, double EE, double PARA6) {
	Ea2 = EE * PARA6 / (fabs(epsun / epsu) + PARA6);
}
void FRCC::Dfunc(double r, double x, double n) {
	if (r != 1) {
		D = 1.0 + (n - r / (r - 1.0)) * (x)+pow(x, r) / (r - 1.0);
	}
	else {
		D = 1.0 + (n - 1.0 + log(x)) * (x);
	}
}
void FRCC::TstressAndTtangent(double epsI, double fI, double EI, double epsF, double fF, double EF, double eps) {

	double RR, AA, EsecAR;
	EsecAR = (fF - fI) / (epsF - epsI);
	if (fabs(EsecAR) > DBL_MAX || fabs(EsecAR) < DBL_MIN) {
		EsecAR = EI;
		RR = 0.0;
		AA = 0.0;
	}
	else {
		RR = (EF - EsecAR) / (EsecAR - EI);
		AA = (EsecAR - EI) / pow(fabs(epsF - epsI), RR);
		
	}
	if (RR < 0.0 || RR>100.0 || fabs(AA) > DBL_MAX || fabs(AA) < DBL_MIN || pow(fabs(epsF - epsI), RR) == 0.0 || pow(fabs(eps - epsI), RR) > DBL_MAX || pow(fabs(eps - epsI), RR) < DBL_MIN) {
		Tstress = fI + (eps - epsI) * EsecAR;
		Ttangent = EsecAR;
	}
	else {
		
		Tstress = fI + (eps - epsI) * (EI + AA * pow((fabs(eps - epsI)), RR));
		Ttangent = EI + AA * (RR + 1) * pow((fabs(eps - epsI)), RR);
	}
}
void FRCC::fnewEnew(double epsMax, double epsul, double sigMax, double mlsig) {

	fnew = mlsig * sigMax;
	Enew = fnew / (epsMax - epsul);
}
void FRCC::EnvComp(double eps, double epsc0) {
	if (eps >= ecu - epsc0) {
		Dfunc(rc, fabs((eps - epsc0) / ecp), nc);

		Tstress = fcp * nc * fabs((eps - epsc0) / ecp) / D;
		Ttangent = Ec * (1.0 - pow(fabs((eps - epsc0) / ecp), rc)) / pow(D, 2.0);
	}
	else {
		Dfunc(rc, fabs((ecu - epsc0) / ecp), nc);

		Ycr = nc * fabs((ecu - epsc0) / ecp) / D;
		Zcr = (1.0 - pow(fabs((ecu - epsc0) / ecp), rc)) / pow(D, 2.0);

		Tstress = fcp * (Ycr + nc * Zcr * fabs((eps - ecu) / ecp));
		Ttangent = Ec * Zcr;
	}
}
void FRCC::EnvTen(double eps, double epst0) {

	if (eps <= et1 - epst0) {
		TstressAndTtangent(-epst0, 0.0, Ec, (et1 - epst0), ft1, Eh, eps);

	}
	else if (eps <= et2 - epst0) {
		Tstress = ft1 + Eh * (eps - (et1 - epst0));
		Ttangent = Eh;
	}
	else if (eps <= etu - epst0) {
		if (fabs(EsecDesce) > DBL_MAX || fabs(EsecDesce) < DBL_MIN) {
			EsecDesce = Ec;
			RRR = 0.0;
			AAA = 0.0;
		}
		if (RRR<0.0 || RRR>100.0 || fabs(AAA)>DBL_MAX || fabs(AAA) < DBL_MIN || pow(fabs(etu - et2), RRR)==0.0 || pow(fabs(etu - et2), RRR)>DBL_MAX || pow(fabs(etu - et2), RRR)<DBL_MIN) {
			Tstress = ft2 + (eps - (et2 - epst0)) * EsecDesce;
			Ttangent = EsecDesce;
		}
		else {
			Tstress = ft2 + (eps - (et2 - epst0)) * (Ec + AAA * pow((fabs(eps - (et2 - epst0))), RRR));
			Ttangent = Ec + AAA * (RRR + 1) * pow((fabs(eps - (et2 - epst0))), RRR);
		}
	}
	else if (eps <= Xcrk) {
		Tstress = -stu * (eps - (etu - epst0)) + ftu;
		Ttangent = -stu;
	}
	else {
		Tstress = 0.0;
		Ttangent = 0.0;
	}
}

double FRCC::getStress()
{
	return Tstress;
}

double FRCC::getStrain()
{
	return Tstrain;
}

double FRCC::getTangent()
{
	return Ttangent;
}

int FRCC::commitState()
{
	Cepstul = Tepstul;
	CsigtMax = TsigtMax;
	Cepscul = Tepscul;
	CsigcMax = TsigcMax;
	CepstMax = TepstMax;
	CepscMax = TepscMax;
	CepstMaxStar = TepstMaxStar;
	CepscMaxStar = TepscMaxStar;
	CsigtMaxStar = TsigtMaxStar;
	CsigcMaxStar = TsigcMaxStar;
	CepstulStar = TepstulStar;
	CepsculStar = TepsculStar;
	CsigtulStar = TsigtulStar;
	CsigculStar = TsigculStar;
	CepsculStar1 = TepsculStar1;
	CsigculStar1 = TsigculStar1;
	CepsculStar2 = TepsculStar2;
	CsigculStar2 = TsigculStar2;
	CepsculStar3 = TepsculStar3;
	CsigculStar3 = TsigculStar3;
	CepstulStar1 = TepstulStar1;
	CsigtulStar1 = TsigtulStar1;
	CepstulStar2 = TepstulStar2;
	CsigtulStar2 = TsigtulStar2;
	CepstulStar3 = TepstulStar3;
	CsigtulStar3 = TsigtulStar3;
	CepstRe = TepstRe;
	CepscRe = TepscRe;
	Cepst0 = Tepst0;
	Cepsc0 = Tepsc0;
	CTI = TTI;
	CFlagT = TFlagT;
	CFlagC = TFlagC;
	Cid1 = Tid1;
	Cid2 = Tid2;
	Cid3 = Tid3;
	Cid4 = Tid4;
	Cid2_1 = Tid2_1;
	Cid4_1 = Tid4_1;
	// State variables
	Cstrain = Tstrain;
	Cstress = Tstress;
	Ctangent = Ttangent;
	return 0;
}

int FRCC::revertToLastCommit()
{
	Tepstul = Cepstul;
	TsigtMax = CsigtMax;
	Tepscul = Cepscul;
	TsigcMax = CsigcMax;
	TepstMax = CepstMax;
	TepscMax = CepscMax;
	TepstMaxStar = CepstMaxStar;
	TepscMaxStar = CepscMaxStar;
	TsigtMaxStar = CsigtMaxStar;
	TsigcMaxStar = CsigcMaxStar;
	TepstulStar = CepstulStar;
	TepsculStar = CepsculStar;
	TsigtulStar = CsigtulStar;
	TsigculStar = CsigculStar;
	TepsculStar1 = CepsculStar1;
	TsigculStar1 = CsigculStar1;
	TepsculStar2 = CepsculStar2;
	TsigculStar2 = CsigculStar2;
	TepsculStar3 = CepsculStar3;
	TsigculStar3 = CsigculStar3;
	TepstulStar1 = CepstulStar1;
	TsigtulStar1 = CsigtulStar1;
	TepstulStar2 = CepstulStar2;
	TsigtulStar2 = CsigtulStar2;
	TepstulStar3 = CepstulStar3;
	TsigtulStar3 = CsigtulStar3;
	TepstRe = CepstRe;
	TepscRe = CepscRe;
	Tepst0 = Cepst0;
	Tepsc0 = Cepsc0;
	TTI = CTI;
	TFlagT = CFlagT;
	TFlagC = CFlagC;
	Tid1 = Cid1;
	Tid2 = Cid2;
	Tid3 = Cid3;
	Tid4 = Cid4;
	Tid2_1 = Cid2_1;
	Tid4_1 = Cid4_1;
	// Recompute trial stress and tangent
	Tstrain = Cstrain;
	Tstress = Cstress;
	Ttangent = Ctangent;

	return 0;
}

int FRCC::revertToStart()
{
	// History variables
	Cepstul = 0.0;
	CsigtMax = 0.0;
	Cepscul = 0.0;
	CsigcMax = 0.0;
	CepstMax = 0.0;
	CepscMax = 0.0;
	CepstMaxStar = 0.0;
	CepscMaxStar = 0.0;
	CsigtMaxStar = 0.0;
	CsigcMaxStar = 0.0;
	CepstulStar = 0.0;
	CepsculStar = 0.0;
	CsigtulStar = 0.0;
	CsigculStar = 0.0;
	CepsculStar1 = 0.0;
	CsigculStar1 = 0.0;
	CepsculStar2 = 0.0;
	CsigculStar2 = 0.0;
	CepsculStar3 = 0.0;
	CsigculStar3 = 0.0;
	CepstulStar1 = 0.0;
	CsigtulStar1 = 0.0;
	CepstulStar2 = 0.0;
	CsigtulStar2 = 0.0;
	CepstulStar3 = 0.0;
	CsigtulStar3 = 0.0;
	CepstRe = 0.0;
	CepscRe = 0.0;
	Cepst0 = 0.0;
	Cepsc0 = 0.0;
	CTI = 0.0;
	CFlagT = 0;
	CFlagC = 0;
	Cid1 = 0;
	Cid2 = 0;
	Cid3 = 0;
	Cid4 = 0;
	Cid2_1 = 0;
	Cid4_1 = 0;
	Tepstul = 0.0;
	TsigtMax = 0.0;
	Tepscul = 0.0;
	TsigcMax = 0.0;
	TepstMax = 0.0;
	TepscMax = 0.0;
	TepstMaxStar = 0.0;
	TepscMaxStar = 0.0;
	TsigtMaxStar = 0.0;
	TsigcMaxStar = 0.0;
	TepstulStar = 0.0;
	TepsculStar = 0.0;
	TsigtulStar = 0.0;
	TsigculStar = 0.0;
	TepsculStar1 = 0.0;
	TsigculStar1 = 0.0;
	TepsculStar2 = 0.0;
	TsigculStar2 = 0.0;
	TepsculStar3 = 0.0;
	TsigculStar3 = 0.0;
	TepstulStar1 = 0.0;
	TsigtulStar1 = 0.0;
	TepstulStar2 = 0.0;
	TsigtulStar2 = 0.0;
	TepstulStar3 = 0.0;
	TsigtulStar3 = 0.0;
	TepstRe = 0.0;
	TepscRe = 0.0;
	Tepst0 = 0.0;
	Tepsc0 = 0.0;
	TTI = 0.0;
	TFlagT = 0;
	TFlagC = 0;
	Tid1 = 0;
	Tid2 = 0;
	Tid3 = 0;
	Tid4 = 0;
	Tid2_1 = 0;
	Tid4_1 = 0;
	// State variables
	Cstrain = 0.0;
	Cstress = 0.0;
	Ctangent = Ec;

	Tstrain = 0.0;
	Tstress = 0.0;
	Ttangent = Ec;

	return 0;
}

UniaxialMaterial* FRCC::getCopy()
{
	FRCC* theCopy = new FRCC(this->getTag(),
		Ec, et1, ft1, et2, ft2, etu, ftu, stu, paraT1, paraT2, paraT3, paraT4, paraT5, paraT6,
		ecp, fcp, ecu, paraC1, paraC2, paraC3, paraC4, paraC5, paraC6, rc);

	// Converged history variables
	theCopy->Cepstul = Cepstul;
	theCopy->CsigtMax = CsigtMax;
	theCopy->Cepscul = Cepscul;
	theCopy->CsigcMax = CsigcMax;
	theCopy->CepstMax = CepstMax;
	theCopy->CepscMax = CepscMax;
	theCopy->CepstMaxStar = CepstMaxStar;
	theCopy->CepscMaxStar = CepscMaxStar;
	theCopy->CsigtMaxStar = CsigtMaxStar;
	theCopy->CsigcMaxStar = CsigcMaxStar;
	theCopy->CepstulStar = CepstulStar;
	theCopy->CepsculStar = CepsculStar;
	theCopy->CsigtulStar = CsigtulStar;
	theCopy->CsigculStar = CsigculStar;
	theCopy->CepsculStar1 = CepsculStar1;
	theCopy->CsigculStar1 = CsigculStar1;
	theCopy->CepsculStar2 = CepsculStar2;
	theCopy->CsigculStar2 = CsigculStar2;
	theCopy->CepsculStar3 = CepsculStar3;
	theCopy->CsigculStar3 = CsigculStar3;
	theCopy->CepstulStar1 = CepstulStar1;
	theCopy->CsigtulStar1 = CsigtulStar1;
	theCopy->CepstulStar2 = CepstulStar2;
	theCopy->CsigtulStar2 = CsigtulStar2;
	theCopy->CepstulStar3 = CepstulStar3;
	theCopy->CsigtulStar3 = CsigtulStar3;
	theCopy->CepstRe = CepstRe;
	theCopy->CepscRe = CepscRe;
	theCopy->Cepst0 = Cepst0;
	theCopy->Cepsc0 = Cepsc0;
	theCopy->CTI = CTI;
	theCopy->CFlagT = CFlagT;
	theCopy->CFlagC = CFlagC;
	theCopy->CFlagT = CFlagT;
	theCopy->CFlagC = CFlagC;
	theCopy->Cid1 = Cid1;
	theCopy->Cid2 = Cid2;
	theCopy->Cid3 = Cid3;
	theCopy->Cid4 = Cid4;
	theCopy->Cid2_1 = Cid2_1;
	theCopy->Cid4_1 = Cid4_1;
	// Converged state variables
	theCopy->Cstrain = Cstrain;
	theCopy->Cstress = Cstress;
	theCopy->Ctangent = Ctangent;

	return theCopy;
}

int FRCC::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;
	static Vector data(67);
	data(0) = this->getTag();

	// Input varibles
	data(1) = Ec;
	data(2) = et1;
	data(3) = ft1;
	data(4) = et2;
	data(5) = ft2;
	data(6) = etu;
	data(7) = ftu;
	data(8) = stu;
	data(9) = paraT1;
	data(10) = paraT2;
	data(11) = paraT3;
	data(12) = paraT4;
	data(13) = paraT5;
	data(14) = paraT6;
	data(15) = ecp;
	data(16) = fcp;
	data(17) = ecu;
	data(18) = paraC1;
	data(19) = paraC2;
	data(20) = paraC3;
	data(21) = paraC4;
	data(22) = paraC5;
	data(23) = paraC6;
	data(24) = rc;

	// Converged varibles
	data(25) = Cepstul;
	data(26) = CsigtMax;
	data(27) = Cepscul;
	data(28) = CsigcMax;
	data(29) = CepstMax;
	data(30) = CepscMax;
	data(31) = CepstMaxStar;
	data(32) = CepscMaxStar;
	data(33) = CsigtMaxStar;
	data(34) = CsigcMaxStar;
	data(35) = CepstulStar;
	data(36) = CepsculStar;
	data(37) = CsigtulStar;
	data(38) = CsigculStar;
	data(39) = CepsculStar1;
	data(40) = CsigculStar1;
	data(41) = CepsculStar2;
	data(42) = CsigculStar2;
	data(43) = CepsculStar3;
	data(44) = CsigculStar3;
	data(45) = CepstulStar1;
	data(46) = CsigtulStar1;
	data(47) = CepstulStar2;
	data(48) = CsigtulStar2;
	data(49) = CepstulStar3;
	data(50) = CsigtulStar3;
	data(51) = CepstRe;
	data(52) = CepscRe;
	data(53) = Cepst0;
	data(54) = Cepsc0;
	data(55) = CTI;
	data(56) = CFlagT;
	data(57) = CFlagC;
	data(58) = Cid1;
	data(59) = Cid2;
	data(60) = Cid3;
	data(61) = Cid4;
	data(62) = Cid2_1;
	data(63) = Cid4_1;
	data(64) = Cstrain;
	data(65) = Cstress;
	data(66) = Ctangent;

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0)
		opserr << "FRCC::sendSelf() - failed to send data\n";

	return res;
}

int FRCC::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int res = 0;
	static Vector data(67);
	res = theChannel.recvVector(this->getDbTag(), commitTag, data);

	if (res < 0) {
		opserr << "FRCC::recvSelf() - failed to receive data\n";
		this->setTag(0);
	}
	else {
		this->setTag(int(data(0)));

		// input varibles
		Ec = data(1);
		et1 = data(2);
		ft1 = data(3);
		et2 = data(4);
		ft2 = data(5);
		etu = data(6);
		ftu = data(7);
		stu = data(8);
		paraT1 = data(9);
		paraT2 = data(10);
		paraT3 = data(11);
		paraT4 = data(12);
		paraT5 = data(13);
		paraT6 = data(14);
		ecp = data(15);
		fcp = data(16);
		ecu = data(17);
		paraC1 = data(18);
		paraC2 = data(19);
		paraC3 = data(20);
		paraC4 = data(21);
		paraC5 = data(22);
		paraC6 = data(23);
		rc = data(24);

		// Converged varibles
		Cepstul = data(25);
		CsigtMax = data(26);
		Cepscul = data(27);
		CsigcMax = data(28);
		CepstMax = data(29);
		CepscMax = data(30);
		CepstMaxStar = data(31);
		CepscMaxStar = data(32);
		CsigtMaxStar = data(33);
		CsigcMaxStar = data(34);
		CepstulStar = data(35);
		CepsculStar = data(36);
		CsigtulStar = data(37);
		CsigculStar = data(38);
		CepsculStar1 = data(39);
		CsigculStar1 = data(40);
		CepsculStar2 = data(41);
		CsigculStar2 = data(42);
		CepsculStar3 = data(43);
		CsigculStar3 = data(44);
		CepstulStar1 = data(45);
		CsigtulStar1 = data(46);
		CepstulStar2 = data(47);
		CsigtulStar2 = data(48);
		CepstulStar3 = data(49);
		CsigtulStar3 = data(50);
		CepstRe = data(51);
		CepscRe = data(52);
		Cepst0 = data(53);
		Cepsc0 = data(54);
		CTI = data(55);
		CFlagT = data(56);
		CFlagC = data(57);
		Cid1 = data(58);
		Cid2 = data(59);
		Cid3 = data(60);
		Cid4 = data(61);
		Cid2_1 = data(62);
		Cid4_1 = data(63);
		Cstrain = data(64);
		Cstress = data(65);
		Ctangent = data(66);

		// State variables from last converged state
		Tstrain = Cstrain;
		Tstress = Cstress;
		Ttangent = Ctangent;
	}
	return res;
}

void FRCC::Print(OPS_Stream& s, int flag)
{
	s << "FRCC, tag: " << this->getTag() << endln;
	s << "  Ec:" << Ec << endln;
	s << "  et1:" << et1 << endln;
	s << "  ft1:" << ft1 << endln;
	s << "  et2:" << et2 << endln;
	s << "  ft2:" << ft2 << endln;
	s << "  etu:" << etu << endln;
	s << "  ftu:" << ftu << endln;
	s << "  stu:" << stu << endln;
	s << "  paraT1:" << paraT1 << endln;
	s << "  paraT2:" << paraT2 << endln;
	s << "  paraT3:" << paraT3 << endln;
	s << "  paraT4:" << paraT4 << endln;
	s << "  paraT5:" << paraT5 << endln;
	s << "  paraT6:" << paraT6 << endln;
	s << "  ecp: " << ecp << endln;
	s << "  fcp:" << fcp << endln;
	s << "  ecu:" << ecu << endln;
	s << "  paraC1:" << paraC1 << endln;
	s << "  paraC2:" << paraC2 << endln;
	s << "  paraC3:" << paraC3 << endln;
	s << "  paraC4:" << paraC4 << endln;
	s << "  paraC5:" << paraC5 << endln;
	s << "  paraC6:" << paraC6 << endln;
	s << "  rc:" << rc << endln;
}