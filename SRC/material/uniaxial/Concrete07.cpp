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
                                                                        
// $Revision: 1.3 $
// $Date: 2007-06-28 22:16:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete07.cpp,v $

// Written: Jon Waugh, Iowa State University
// Created: 10/2006
//
// Description: This file contains the class definition for Uniaxial material Concrete07
//				A implementation of the Chang & Mander Concrrete model from 1994.

#include <Concrete07.h>
#include <UniaxialMaterial.h>
#include <math.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <float.h>
#include <iostream>
#include <elementAPI.h>

void* OPS_Concrete07()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 9) {
	opserr << "WARNING: Insufficient arguments\n";
	opserr << "Want: uniaxialMaterial Concrete07 tag? ";
	opserr << "fpc? epsc0? Ec? fpt? epst0? xcrp? xcrn? r?\n";
	return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING invalid tag\n";
	return 0;
    }

    double data[8];
    numdata = 8;
    if (OPS_GetDoubleInput(&numdata,data)) {
	opserr << "WARNING invalid double data\n";
	return 0;
    }

    UniaxialMaterial* mat = new Concrete07(tag,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7]);
    if (mat == 0) {
	opserr << "WARNING: failed to create Concrete07 material\n";
	return 0;
    }

    return mat;
}

Concrete07::Concrete07 (int tag, double FPC, double EPSC0, double EC, double FPT, double EPST0, double XCRP, double XCRN, double R)
:UniaxialMaterial(tag, MAT_TAG_Concrete07), fpc(FPC), epsc0(EPSC0), Ec(EC), fpt(FPT), epst0(EPST0), xcrp(XCRP), xcrn(XCRN), r(R) {

	// Calculate the variables that are needed to define the envelopw
	nn = (Ec*epsc0)/fpc;
	np = (Ec*epst0)/fpt;

	double y(0),z(0);
	
	calculateYandZ(xcrn,y,z,nn);

	xsp = xcrn - y/(nn*z);

	calculateYandZ(xcrp,y,z,np);

	xcrk = xcrp - y/(np*z);

	e0 = 0;

	// Set all history and state variables to initial values
	this->revertToStart();
}

Concrete07::Concrete07() : UniaxialMaterial(0,MAT_TAG_Concrete07)
{
	opserr << "WARNING: Reguire input of tag, fpc, epsc0, Ec, fpt, epst0, xcrp, xcrn, deltaFcu, r\n";
}

Concrete07::~Concrete07()
{
	// No dynamic variables are used so a destructor is not required.
}


int Concrete07::setTrialStrain(double strain, double strainRate)
{
	// Reset History variables to last converged state
	TminStrain = CminStrain;
	TminStress = CminStress;
	TUnloadStrain = CUnloadStrain;
	TUnloadStress = CUnloadStress;
	TUnloadStiffness = CUnloadStiffness;
	TmaxStrain = CmaxStrain;
	TmaxStress = CmaxStress;
	TReloadStrain = CReloadStrain;
	TReloadStress = CReloadStress;
	Tloading = Cloading;
	Tcracked = Ccracked;
	Trule = Crule;

	// Set trial strain
	Tstrain = strain;
	
	// Determine change in strain
	double dStrain = strain - Cstrain;

	// Calculate the trial state given the trial strain
	determineTrialState(dStrain);

	return 0;
}

int Concrete07::setTrial(double strain, double &stress, double &tangent, double strainRate)
{
	// Reset History variables to last converged state
	TminStrain = CminStrain;
	TminStress = CminStress;
	TUnloadStrain = CUnloadStrain;
	TUnloadStress = CUnloadStress;
	TUnloadStiffness = CUnloadStiffness;
	TmaxStrain = CmaxStrain;
	TmaxStress = CmaxStress;
	TReloadStrain = CReloadStrain;
	TReloadStress = CReloadStress;
	T13Zero = C13Zero;
	Tloading = Cloading;
	Tcracked = Ccracked;
	Trule = Crule;

	// Set trial strain
	Tstrain = strain;

	// Determine change in strain
	double dStrain = strain - Cstrain;

	// Calculate the trial state given the trial strain
	determineTrialState(dStrain);

	stress = Tstress;
	tangent = Ttangent;

	return 0;
}

void Concrete07::calculateYandZ(double x, double& y, double& z, double n)
{
	double D;
	if (r == 1)
		D = 1+(n-1+log(x))*x;
	else
		D = 1+(n-r/(r-1))*x+pow(x,r)/(r-1);

	y = n*x/D;
	z = (1-pow(x,r))/(pow(D,2));

	return;
}


void Concrete07::calculateStressTransition(double& fc, double& Et, double ec, double eI, double fI, double EI, double eF, double fF, double EF, int rule)
{
	double er, ea, eb, fa, fb;
	int dir; 

	// where er is the strain at the breakpoint.  This is the intersetcion of the two straight 
	// lines defined by the beginning and ending points and slopes.

	// Calculate the strain at the breakpoint between the two lines
	er = (EI*eI - EF*eF - fI + fF)/(EI-EF);
	ea = (er + eI) /2;
	eb = (eF + er) /2;
	fa = EI*(ea - eI)+fI;
	fb = EF*(eb - eF)+fF;

	// Determine if the strain is above or below the breakpoint and calculate the stress and stiffness accordingly

	if (eI < eF)
	{
		dir = 1;

		// Check to ensure that point R is between our two strains.
		if (er >= eF)
		{
			Et = fabs((fF-fI)/(eF-eI));
			fc = Et*(ec-eI)+fI;
			
			return;
		}
	}

	else
	{
		dir = 2;
		// Check to ensure that point R is between our two strains.
		if (eF >= er)
		{
			Et = fabs((fF-fI)/(eF-eI));
			fc = Et*(ec-eI)+fI;

			return;
		}
	}

	switch (dir)
	{
		case 1:
			if (ec <= ea)
			{
				fc = EI*(ec - eI) + fI;
				Et = EI;
			}

			else if (ec <= eb)
			{
				Et = (fb - fa)/(eb-ea);
				fc = Et*(ec - ea) + fa;
			}

			else
			{
				fc = EF*(ec - eF) + fF;
				Et = EF;
			}

			break;

		case 2:
			if (ec >= ea)
			{
				fc = EI*(ec - eI) + fI;
				Et = EI;
			}

			else if (ec >= eb)
			{
				Et = (fb - fa)/(eb-ea);
				fc = Et*(ec - ea) + fa;
			}

			else
			{
				fc = EF*(ec - eF) + fF;
				Et = EF;
			}

			break;
	}

	return;
}


void Concrete07::calculate13Stress(double &fc, double &Et, double ec, double eI, double eF, double fF, double EF)
{
	double A;				// Equation parameter
	double R(0);			// Equation Parameter
	double ESEC;			// Secant Modulus
	double fI = 0;			// Initial Stress
	double EI = 0;			// Youngs Modulus

	ESEC = (fF-fI)/(eF-eI);

	if ( EI/ESEC >= 0.985 && EI/ESEC < 1.015)
		R = 0;
	else
		R = fabs((EF-ESEC)/(ESEC-EI));
	 if (R > 100)
	{
		calculateStressTransition(fc, Et, ec, eI, 0.0, 0.25*ESEC, eF, fF, EF, 666);
		
		return;
	}

	if (eF/eI > 0.9999 && eF/eI < 1.0001)
		R = 0;
	
	A = (ESEC-EI)/pow(fabs(eF-eI),R);

	if (A > pow(10.0,300.0))
	{
		calculateStressTransition(fc, Et, ec, eI, 0.0, 0.25*ESEC, eF, fF, EF, 666);

		return;
	}


	fc = fI+(ec - eI)*(EI+A*(pow(fabs(ec-eI),R)));
	Et = EI+A*(R+1)*pow(fabs(ec-eI),R);

	return;
}
void Concrete07::envelope(double x, double& fc, double& Et, int flag)
{
	double y,z;
	
	if (flag >= 0)
	{
			
		if (x<xcrp)
		{
			calculateYandZ(x,y,z,np);

			fc = fpt*y;
			Et = Ec*z;
			Trule = 2;
		}

		else if (x<=xcrk)
		{
			calculateYandZ(xcrp,y,z,np);

			fc = fpt*(y+np*z*(x-xcrp));
			Et = Ec*z;
			Trule = 2;
		}

		else
		{
			fc = 0.0;
			Et = 0.0;
			Trule = 6;
		}
	}

	else
	{
		if (x<xcrn)
		{
			calculateYandZ(x,y,z,nn);

			fc = fpc*y;
			Et = Ec*z;
			Trule = 1;
		}

		else if (x<=xsp)
		{
			calculateYandZ(xcrn,y,z,nn);

			fc = fpc*(y+nn*z*(x-xcrn));
			Et = Ec*z;
			Trule = 1;
		}

		else 
		{
			fc = 0.0;
			Et = 0.0;
			Trule = 5;
		}
	}

	return;
}

void Concrete07::determineTrialState(double dStrain)
{
	double eunn;			// unloading strain from the compression envelope
	double funn;			// unloading stress from the compression envelope 
	double Esecn;			// secant modulus at loading point in compression
	double epln;			// plastic strain in compression
	double eunp;			// reloading strain from the tension envelope
	double funp;			// reloading stress from the tension envelope	
	double Esecp;			// secant modulus at loading point in tension
	double eplp;			// plastic strain in tension
	double x;

	// Calculate all the variables
	eunn = TminStrain;
	funn = TminStress;

	Esecn = Ec*((funn/(Ec*epsc0)+0.57)/(eunn/epsc0+0.57));
	epln = eunn-funn/Esecn;

	eunp = TmaxStrain;
	e0 = 0.0;	

	//if (fabs((eunp-e0)/epst0)<(eunn/epsc0))
	if (eunp == 0 && eunn/epsc0 <= xsp)
	{
		int ruleStore = Trule;
		double Edum;
		eunp = eunn/epsc0*epst0;
		envelope(eunp/epst0, funp, Edum, 1);
		Trule = ruleStore;
		TmaxStrain = eunp;
		TmaxStress = funp;
	}
	else if (eunn/epsc0 > xsp)
		Tcracked = true;
	else
	{
		funp = TmaxStress;
	}

	Esecp = Ec*((funp/(Ec*epst0)+0.67)/((eunp-e0)/epst0+0.67));
	eplp = eunp - funp/Esecp;
	if (eplp < 0)
		eplp = 0;

	// Calculate x for the given strain;
	
	if (Tstrain >= 0.0)
	{
		x = fabs((Tstrain-e0)/epst0);
		if (x > xcrk)
			Tcracked = true;
	}

	else
	{
		x = Tstrain/epsc0;
	}

	if (fabs(dStrain) <= DBL_EPSILON)    // ignore trivial strain change
	{
		Tstress = Cstress;
		return;
	}

	if (Tloading == 0)
	{
		if (Tstrain >= 0)
		{
			
			envelope(x, Tstress, Ttangent, 1);

			Tloading = -1;

			TmaxStrain = Tstrain;
			TmaxStress = Tstress;

		}

		else
		{
			envelope(x, Tstress, Ttangent, -1);

			Tloading = 1;

			TminStrain = Tstrain;
			TminStress = Tstress;

		}

		return;
	}

	if (Tloading > 0)				// Previously loading the concrete fibers
	{
		if (dStrain < 0.0)			// Continue loading the concrete (compression)
		{
			if (Tstrain <= eunn)
			{
				
				envelope(x, Tstress, Ttangent, -1);

				TminStrain = Tstrain;
				TminStress = Tstress;

				return;
			}
				
			else if (Tstrain < eplp && !Tcracked)
			{
				if (Trule == 71)				// we are on the transition curve for reloading from a partial unloading in compression
				{
					double eron;			// strain at which reversal occurred in a partial unloading
					double fron;			// stress at which reversal occurred in a partial unloading

					eron = TReloadStrain;
					fron = TReloadStress;

					Ttangent = (funn - fron)/(eunn - eron);
					Tstress = Ttangent * (Tstrain - eron) + fron;

					Trule = 71;

					return;
				}

				else if (Trule == 11)				// we are on the transition curve for parital reloading after a partial unloading
				{
					double er = TReloadStrain;	// strain when reversal occurred.
					double fr = TReloadStress;	// stress when reversal occurred.
					double Eplp;				// Modulus at plastic point in tension
					double fb;					// target stress on rule 10
					double eb;					// target strain on rule 10;
					double Eb;					// slope at target point (eb,fb)
					double Enewn;				// reloading stiffness in compression;

					Eplp = Ec/(pow(fabs((eunp-e0)/epst0),1.1)+1.00);
					Enewn = funn/(eunn-epln);;
					eb = eunn - (er-epln)/(eunp-epln)*(eunn-eplp);

					if (Tstrain > eb)			// are we still on the transition curve
					{

						calculateStressTransition(fb, Eb, eb, eplp, 0.0, Eplp, eunn, funn, Enewn, 10);

						Ttangent = (fb - fr)/(eb - er);
						Tstress = Ttangent*(Tstrain - er) + fr;

						Trule = 11;

						return;
					}

					else						// we are back on the connecting curve;
					{
						calculateStressTransition(Tstress, Ttangent, Tstrain, eplp, 0.0, Eplp, eunn, funn, Enewn, 10);

						Trule = 10;

						return;
					}
				}


				else 							// we are reloading from a complete unloading
				{
					double Enewn;			// reloading stiffness in compression;
					double Eplp;			// Modulus at plastic point in tension

					if (eunn == 0 && funn == 0)
					{
						double Edum;
						eunn = -0.00002;
						envelope(eunn/epsc0, funn, Edum, -1);
						TminStrain = eunn;
						TminStress = funn;
						Esecn = Ec*((funn/(Ec*epsc0)+0.57)/(eunn/epsc0+0.57));
						epln = eunn-funn/Esecn;
					}

					Enewn = funn/(eunn-epln);
					Eplp = Ec/(pow(fabs((eunp-e0)/epst0),1.1)+1.00);
					calculateStressTransition(Tstress, Ttangent, Tstrain, eplp, 0.0, Eplp, eunn, funn, Enewn, 10);

					Trule = 10;

					return;
				}
			}

			else if (!Tcracked)			// We are reloading in the tension side
			{
				if (Trule == 4)
				{
					double Eplp;			// Modulus at plastic point in tension

					Eplp = Ec/(pow(fabs((eunp-e0)/epst0),1.1)+1.00);

					calculateStressTransition(Tstress, Ttangent, Tstrain, eunp, funp, Ec, eplp, 0.0, Eplp, 4);

					Trule = 4;

					return;
				}

				else if (Trule == 811 && Tstrain > TUnloadStrain)
				{
					double Enewps;			// stiffness at strain reloading from tension envelope after unloading
					double erop;			// strain at last reversal;
					double frop;			// stress at last reversal;

					erop = TUnloadStrain;
					frop = TUnloadStress;
					Enewps = (funp - frop)/(eunp - erop);

					Ttangent = Enewps;
					Tstress = Ttangent*(Tstrain - erop) + frop;

					Trule = 811;

					return;
				}


				else if (Trule == 11)
				{
					double er = TReloadStrain;	// strain when reversal occurred.
					double fr = TReloadStress;	// stress when reversal occurred.
					double Eplp;				// Modulus at plastic point in tension
					double fb;					// target stress on rule 10
					double eb;					// target strain on rule 10;
					double Eb;					// slope at target point (eb,fb)
					double Enewn;				// reloading stiffness in compression;

					Eplp = Ec/(pow(fabs((eunp-e0)/epst0),1.1)+1.00);
					Enewn = funn/(eunn-epln);;
					eb = eunn - (er-epln)/(eunp-epln)*(eunn-eplp);

					if (Tstrain > eb)			// are we still on the transition curve
					{

						calculateStressTransition(fb, Eb, eb, eplp, 0.0, Eplp, eunn, funn, Enewn, 10);

						Ttangent = (fb - fr)/(eb - er);
						Tstress = Ttangent*(Tstrain - er) + fr;

						Trule = 11;

						return;
					}

					else						// we are back on the connecting curve;
					{
						calculateStressTransition(Tstress, Ttangent, Tstrain, eplp, 0.0, Eplp, eunn, funn, Enewn, 10);

						Trule = 10;

						return;
					}
				}
				else 
				{
					double Eplp;			// Modulus at plastic point in tension

					Eplp = Ec/(pow(fabs((eunp-e0)/epst0),1.1)+1.00);

					calculateStressTransition(Tstress, Ttangent, Tstrain, eunp, funp, Ec, eplp, 0.0, Eplp, 4);

					Trule = 4;

					return;
				}
			}

			else if (Trule == 71)				// we are on the transition curve for reloading from a partial unloading in compression
			{
				double eron;			// strain at which reversal occurred in a partial unloading
				double fron;			// stress at which reversal occurred in a partial unloading
				double Enewns;			// stiffness at strain unloading from compression envelope after reloading

				eron = TReloadStrain;
				fron = TReloadStress;
				Enewns = (funn-fron)/(eunn-eron);

				Ttangent = Enewns;
				Tstress = Ttangent * (Tstrain - eron) + fron;

				Trule = 71;

				return;
			}

			else if (Trule == 15)
			{
				double ea = T13Strain;			// Strain at unloading on rule 13
				double fa = T13Stress;			// Stress at unloading on rule 13
				double er = TReloadStrain;		// Strain at last reloading
				double fr = TReloadStress;		// Stress at last reloading

				if (Tstrain > ea)
				{
					Ttangent = (fa - fr) / (ea - er);
					Tstress = Ttangent * (Tstrain - er) + fr;

					Trule = 15;

					return;
				}

				else
				{
					double Enewn;					// reloading stiffness in compression;
					double Edum;

					Enewn = funn/(eunn-epln);
					
					Edum = (funn) / (eunn - T13Zero);
					calculate13Stress(Tstress, Ttangent, Tstrain, T13Zero, eunn, funn, Enewn);
					
					Trule = 13;

					return;
				}
			}

			else						// we are reloading after cracking, require gradual crack closure
			{
				if (eunn/epsc0 >= xsp)			// See if we have spalled and make a quick return
				{
					Tstress = 0.0;
					Ttangent = 0.0;
					Trule = 5;
					Tcracked = true;
					
					return;
				}

				double er;				// strain at which reloading begins after cracking
				double Enewn;			// reloading stiffness in compression;
				double Edum;

				er = T13Zero;
				Enewn = funn/(eunn-epln);

				Edum = (funn) / (eunn - er);

				calculate13Stress(Tstress, Ttangent, Tstrain, er, eunn, funn, Enewn);

				Trule = 13;

				return;
			}
		}

		else if (dStrain > 0.0)
		{
			// Previously Loading, now unloading, strain has reversed.  Need to determine what rule we are on.
			
			if (Trule == 1)
			{
				// Reversing from the compression envelop envelope.  We use rule 3
				Tloading = -1;			// We are now unloading the concrete
				TminStress = Cstress;
				TminStrain = Cstrain;
				eunn = Cstrain;
				funn = Cstress;
				
				double Epln;			// Modulus at plastic point in compression
				double Enewp;			// reloading stiffness in tension;

				Epln = 0.1*Ec*exp(-2*(eunn/epsc0));
				Enewp = funp/(eunp-eplp);
				
				if (Tstrain <= epln)
				{
					calculateStressTransition(Tstress, Ttangent, Tstrain, eunn, funn, Ec, epln, 0.0, Epln, 3);

					Trule = 3;

					return;
				}
				
				else if (Tstrain < eunp && !Tcracked)
				{
					calculateStressTransition(Tstress, Ttangent, Tstrain, epln, 0.0, Epln, eunp, funp, Enewp, 9);

					Trule = 9;

					return;
				}	

				else if (!Tcracked)
				{
					envelope (x, Tstress, Ttangent, 1);
					
					TmaxStrain = Tstrain;
					TmaxStress = Tstress;

					return;
				}

				else			// The concrete has cracked
				{
					Tstress = 0.0;
					Tstrain = 0.0;

					Trule = 6;

					return;
				}
			}

			else if (Trule == 4)
			{
				// Partial reloading in tension zone
				Tloading = -1;

				double erop = Cstrain;
				double frop = Cstress;
				
				TUnloadStrain = erop;
				TUnloadStress = frop;
								
				if (Tstrain < eunp)			// Are in the transition curve for unloading from partial reloading
				{
					Ttangent = (funp - frop)/(eunp - erop);
					Tstress = Ttangent*(Tstrain - erop)+frop;

					Trule = 81;

					return;
				}

				else							// Back on the tension envelope
				{
					envelope(x, Tstress, Ttangent, 1);

					TmaxStrain = Tstrain;
					TmaxStress = Tstress;

					return;
				}
			}

			else if (Trule == 811)
			{
				// Partial reloading in tension zone
				Tloading = -1;

				double erop = TUnloadStrain;
				double frop = TUnloadStress;				

				if (Tstrain < eunp)			// Are in the transition curve for unloading from partial reloading
				{
					Ttangent = (funp - frop)/(eunp - erop);
					Tstress = Ttangent*(Tstrain - erop)+frop;

					Trule = 81;

					return;
				}

				else							// Back on the tension envelope
				{
					envelope(x, Tstress, Ttangent, 1);

					TmaxStrain = Tstrain;
					TmaxStress = Tstress;

					return;
				}
			}

			else if (Trule == 71)
			{
				Tloading = -1;

				double er = Cstrain;	// strain when reversal occurred.
				double fr = Cstress;	// stress when reversal occurred.
				double eron;			// strain at which reversal occurred in a partial unloading
				double fron;			// stress at which reversal occurred in a partial unloading

				eron = TReloadStrain;
				fron = TReloadStress;

				if (Tstrain < eron)
				{
					Ttangent = (funn - fron)/(eunn - eron);
					Tstress = Ttangent * (Tstrain - eron) + fron;

					Trule = 711;

					return;
				}


				double Epln;			// Modulus at plastic point in compression
				double Enewp;			// reloading stiffness in tension;

				Epln = 0.1*Ec*exp(-2*(eunn/epsc0));
				Enewp = funp/(eunp-eplp);
				
				if (Tstrain <= epln)
				{
					calculateStressTransition(Tstress, Ttangent, Tstrain, eunn, funn, Ec, epln, 0.0, Epln, 3);

					Trule = 3;

					return;
				}

				else if (Tstrain < eunp && !Tcracked)
				{
					calculateStressTransition(Tstress, Ttangent, Tstrain, epln, 0.0, Epln, eunp, funp, Enewp, 9);

					Trule = 9;

					return;
				}	

				else if (!Tcracked)
				{
					envelope (x, Tstress, Ttangent, 1);
					
					TmaxStrain = Tstrain;
					TmaxStress = Tstress;

					return;
				}

				else			// The concrete has cracked
				{
					Tstress = 0.0;
					Tstrain = 0.0;

					Trule = 6;

					return;
				}
			}		

			else if (Trule == 10 || Trule == 11)
			{
				// Partial reloading in compression zone
				Tloading = -1;

				double er = Cstrain;	// strain when reversal occurred.
				double fr = Cstress;	// stress when reversal occurred.
				double Eplp;			// Modulus at plastic point in tension
				double Epln;			// Modulus at plastic point in compression
				double fa;				// target stress on rule 9
				double ea;				// target strain on rule 9;
				double Ea;				// slope at target point (ea,fa)
				double Enewp;			// reloading stiffness in tension;

				TUnloadStrain = er;
				TUnloadStress = fr;
				TUnloadStiffness = Ctangent;
				Eplp = Ec/(pow(fabs((eunp-e0)/epst0),1.1)+1.00);
				Epln = 0.1*Ec*exp(-2*(eunn/epsc0));
				Enewp = funp/(eunp-eplp);
				ea = (eunn-er)/(eunn-eplp)*(eunp-epln)+epln;
				
				if (Tstrain < ea)				// On the transition curve 
				{
					calculateStressTransition(fa, Ea, ea, epln, 0.0, Epln, eunp, funp, Enewp, 9);

					//calculateStressTransition(Tstress, Ttangent, Tstrain, er, fr, Ec, ea, fa, Ea);

					Ttangent = (fa - fr)/(ea - er);
					Tstress = Ttangent*(Tstrain - er) + fr;

					Trule = 12;

					return;
				}

				else if (Tstrain < eunp)		// On the connecting curve for unloadinf
				{
					calculateStressTransition(Tstress, Ttangent, Tstrain, epln, 0.0, Epln, eunp, funp, Enewp, 9);

					Trule = 9;

					return;
				}

				else							// On the tension envelope
				{
					envelope(x, Tstress, Ttangent, 1);

					TmaxStrain = Tstrain;
					TmaxStress = Tstress;

					return;
				}

			}

			else							// Unloading after cracking of the concrete has occurred.
			{

				Tloading = -1;

				if (Trule == 5)
				{
					Tstress = 0.0;
					Ttangent = 0.0;
					Tcracked = true;

					Trule = 6;

					return;
				}


				double er = Cstrain;	// strain when reversal occurred.
				double fr = Cstress;	// stress when reversal occurred.
				TUnloadStiffness = Ctangent;
				TUnloadStrain = er;
				TUnloadStress = fr;
				double eb = er - fr/Esecn;

				if (Tstrain < eb)  // Going down to rule 6;
				{
					if (Trule == 13)
					{
						T13Strain = er;
						T13Stress = fr;
					}

					Ttangent = (0.0 - fr)/(eb - er);
					Tstress = Ttangent*(Tstrain - er) + fr;

					Trule = 14;

					return;
				}

				// We are on the ordinate acis, follow rule six

				Tstress = 0.0;
				Ttangent = 0.0;

				Trule = 6;

				return;
			}

		}
	}

	if (Tloading < 0)					// Previously where unloading the concrete
	{
		if (dStrain > 0.0)				// We are continuing to unload the concrete
		{
			if (Tstrain > eunp && !Tcracked)
			{

				envelope(x, Tstress, Ttangent, 1);

				TmaxStrain = Tstrain;
				TmaxStress = Tstress;

				return;
			}

			else if (Tstrain > epln && !Tcracked)
			{
				if (Trule == 81)
				{
					double Enewps;			// stiffness at strain reloading from tension envelope after unloading
					double erop;			// strain at last reversal;
					double frop;			// stress at last reversal;

					erop = TUnloadStrain;
					frop = TUnloadStress;
					Enewps = (funp - frop)/(eunp - erop);

					Ttangent = Enewps;
					Tstress = Ttangent*(Tstrain - erop) + frop;

					Trule = 81;

					return;
				}

				else if (Trule == 12)
				{
					double er = TUnloadStrain;			// strain when unloading began
					double fr = TUnloadStress;			// stress when unloading began
					double Epln;						// Modulus at plastic point from compression side
					double fa;							// target stress on rule 9
					double ea;							// target strain on rule 9
					double Ea;							// slope at targeet point (ea, fa)
					double Enewp;						// reloading stiffness in tension;

					Enewp = funp/(eunp-eplp);
					Epln = 0.1*Ec*exp(-2*(eunn/epsc0));
					ea = (eunn-er)/(eunn-eplp)*(eunp-epln)+epln;

					if (Tstrain < ea)			// we are still on the transition curve
					{
						calculateStressTransition(fa, Ea, ea, epln, 0.0, Epln, eunp, funp, Enewp, 9);

						Ttangent = (fa - fr) / (ea - er);
						Tstress = Ttangent*(Tstrain - er) + fr;
						
						Trule = 12;

						return;
					}

					else						// we are back on the connecting curve
					{
						calculateStressTransition(Tstress, Ttangent, Tstrain, epln, 0.0, Epln, eunp, funp, Enewp, 9);

						Trule = 9;

						return;
					}
				}

				else					// we are unloading from a complete reloading
				{
					double Enewp;						// reloading stiffness in tension;
					double Epln;						// Modulus at plastic point from compression side

					Enewp = funp/(eunp-eplp);
					Epln = 0.1*Ec*exp(-2*(eunn/epsc0));

					calculateStressTransition(Tstress, Ttangent, Tstrain, epln, 0.0, Epln, eunp, funp, Enewp, 9);

					Trule = 9;

					return;
				}
			}

			else if (Tstrain < epln)						// we are unloading off of the compression envelope
			{
				if (Trule == 12)
				{
					double er = TUnloadStrain;			// strain when unloading began
					double fr = TUnloadStress;			// stress when unloading began
					double Epln;						// Modulus at plastic point from compression side
					double fa;							// target stress on rule 9
					double ea;							// target strain on rule 9
					double Ea;							// slope at targeet point (ea, fa)
					double Enewp;						// reloading stiffness in tension;

					Enewp = funp/(eunp-eplp);
					Epln = 0.1*Ec*exp(-2*(eunn/epsc0));
					ea = (eunn-er)/(eunn-eplp)*(eunp-epln)+epln;

					if (Tstrain < ea)			// we are still on the transition curve
					{
						calculateStressTransition(fa, Ea, ea, epln, 0.0, Epln, eunp, funp, Enewp, 9);

						Ttangent = (fa -fr) / (ea - er);
						Tstress = Ttangent*(Tstrain - er) + fr;

						Trule = 12;

						return;
					}

					else						// we are back on the connecting curve
					{
						calculateStressTransition(Tstress, Ttangent, Tstrain, epln, 0.0, Epln, eunp, funp, Enewp, 9);

						Trule = 9;

						return;
					}
				}

				if (Trule == 14)
				{
					double er, fr;
					er = TUnloadStrain;
					fr = TUnloadStress;
					double eb = er - fr/Esecn;

					if (Tstrain < eb)  // Going down to rule 6;
					{
						Ttangent = (0.0 - fr)/(eb - er);
						Tstress = Ttangent*(Tstrain - er) + fr;
						
						Trule = 14;

						return;
					}

					else 
					{
						Tstress = 0.0;
						Ttangent = 0.0;

						Trule = 6;

						return;
					}
				}

				if (Trule == 711 && Tstrain < TReloadStrain)
				{
					double eron;			// strain at which reversal occurred in a partial unloading
					double fron;			// stress at which reversal occurred in a partial unloading

					eron = TReloadStrain;
					fron = TReloadStress;

					Ttangent = (funn - fron)/(eunn - eron);
					Tstress = Ttangent * (Tstrain - eron) + fron;

					Trule = 711;

					return;
				}

				double Epln;			// Modulus at plastic point in compression
				
				Epln = 0.1*Ec*exp(-2*(eunn/epsc0));
				
				calculateStressTransition(Tstress, Ttangent, Tstrain, eunn, funn, Ec, epln, 0.0, Epln, 3);

				Trule = 3;
				
				return;
			}

			else
			{
				if (Trule == 14)
				{
					double er, fr;
					er = TUnloadStrain;
					fr = TUnloadStress;
					double eb = er - fr/Esecn;

					if (Tstrain < eb)  // Going down to rule 6;
					{
						Ttangent = (0.0 - fr)/(eb -er);
						Tstress = Ttangent*(Tstrain - er) + fr;
						Trule = 14;

						return;
					}
				}

				Tstress = 0.0;
				Ttangent = 0.0;

				Trule = 6;

				return;
			}
		}

		if (dStrain < 0.0)				// Previously unloading, now reloading
		{
			// We need to determine what rule we are on

			if (Trule == 2)
			{
				Tloading = 1;

				TmaxStress = Cstress;
				TmaxStrain = Cstrain;
				eunp = Cstrain;
				funp = Cstress;

				double Eplp;			// Modulus at plastic point in tension
				double Enewn;			// reloading stiffness in compression;
				
				if (eunn == 0 && funn == 0)
				{
					double Edum;
					eunn = -0.00002;
					envelope(eunn/epsc0, funn, Edum, -1);
					TminStrain = eunn;
					TminStress = funn;
					Esecn = Ec*((funn/(Ec*epsc0)+0.57)/(eunn/epsc0+0.57));
					epln = eunn-funn/Esecn;
				}
				
				Eplp = Ec/(pow(fabs((eunp-e0)/epst0),1.1)+1.00);
				Enewn = funn/(eunn-epln);

				if (Tstrain >= eplp)
				{
					calculateStressTransition(Tstress, Ttangent, Tstrain, eunp, funp, Ec, eplp, 0.0, Eplp, 4);
					
					Trule = 4;

					return;
				}

				else if (Tstrain >= eunn)
				{
					calculateStressTransition(Tstress, Ttangent, Tstrain, eplp, 0.0, Eplp, eunn, funn, Enewn, 10);

					Trule = 10;

					return;
				}

				else
				{
					envelope(x, Tstress, Ttangent, -1);

					TminStress = Tstress;
					TminStrain = Tstrain;

					return;
				}
			}

			else if (Trule == 3)
			{
				// Partial unloading in compression region

				Tloading = 1;

				double eron = Cstrain;
				double fron = Cstress;
				double Enewns;			// stiffness at strain unloading from compression envelope after reloading

				TReloadStrain = eron;
				TReloadStress = fron;
				Enewns = (funn-fron)/(eunn-eron);

				if (Tstrain > eunn)			// On transition due to partial unloading
				{
					Ttangent = Enewns;
					Tstress = Ttangent*(Tstrain - eron) + fron;

					Trule = 71;

					return;
				}

				else						// Back on the compression envelope
				{
					envelope(x, Tstress, Ttangent, -1);

					TminStress = Tstress;
					TminStrain = Tstrain;

					return;
				}
			}

			else if (Trule == 711)
			{
				Tloading = 1;

				double eron;			// strain at which reversal occurred in a partial unloading
				double fron;			// stress at which reversal occurred in a partial unloading

				eron = TReloadStrain;
				fron = TReloadStress;

				if (Tstrain > eunn)
				{

					Ttangent = (funn - fron)/(eunn - eron);
					Tstress = Ttangent * (Tstrain - eron) + fron;

					Trule = 71;

					return;
				}

				else						// Back on the compression envelope
				{
					envelope(x, Tstress, Ttangent, -1);

					TminStress = Tstress;
					TminStrain = Tstrain;

					return;
				}
			}

			else if (Trule == 81)
			{
				Tloading = 1;

				double Enewps;			// stiffness at strain reloading from tension envelope after unloading
				double erop;			// strain at last reversal;
				double frop;			// stress at last reversal;
				double Eplp;			// Modulus at plastic point in tension
				double Enewn;			// reloading stiffness in compression;

				if (eunn == 0 && funn == 0)
				{
					double Edum;
					eunn = -0.00002;
					envelope(eunn/epsc0, funn, Edum, -1);
					TminStrain = eunn;
					TminStress = funn;
					Esecn = Ec*((funn/(Ec*epsc0)+0.57)/(eunn/epsc0+0.57));
					epln = eunn-funn/Esecn;
				}
				
				Eplp = Ec/(pow(fabs((eunp-e0)/epst0),1.1)+1.00);
				Enewn = funn/(eunn-epln);
				erop = TUnloadStrain;
				frop = TUnloadStress;
				Enewps = (funp - frop)/(eunp - erop);

				if (Tstrain > erop)
				{
					Ttangent = Enewps;
					Tstress = Ttangent*(Tstrain - erop) + frop;

					Trule = 811;

					return;
				}

				else if (Tstrain >= eplp)
				{
					calculateStressTransition(Tstress, Ttangent, Tstrain, eunp, funp, Ec, eplp, 0.0, Eplp, 4);
					
					Trule = 4;

					return;
				}

				else if (Tstrain >= eunn)
				{
					calculateStressTransition(Tstress, Ttangent, Tstrain, eplp, 0.0, Eplp, eunn, funn, Enewn, 10);

					Trule = 10;

					return;
				}

				else
				{
					envelope(x, Tstress, Ttangent, -1);

					TminStress = Tstress;
					TminStrain = Tstrain;

					return;
				}
			}

			else if (Trule == 9 || Trule == 12)
			{

				// Reloading after partial unloading

				Tloading = 1;

				double er = Cstrain;	// strain when reversal occurred.
				double fr = Cstress;	// stress when reversal occurred.
				double Eplp;			// Modulus at plastic point in tension
				double fb;				// target stress on rule 10
				double eb;				// target strain on rule 10;
				double Eb;				// slope at target point (eb,fb)
				double Enewn;			// reloading stiffness in tension;

				TReloadStrain = er;
				TReloadStress = fr;
				Eplp = Ec/(pow(fabs((eunp-e0)/epst0),1.1)+1.00);
				Enewn = funn/(eunn-epln);
				eb = eunn - (er-epln)/(eunp-epln)*(eunn-eplp);

				if (Tstrain > eb)
				{
					calculateStressTransition(fb, Eb, eb, eplp, 0.0, Eplp, eunn, funn, Enewn, 10);

					Ttangent = (fb - fr)/(eb - er);
					Tstress = Ttangent*(Tstrain - er) + fr;

					Trule = 11;

					return;
				}

				else if (Tstrain > eunn)		// On the connecting curve
				{

					calculateStressTransition(Tstress, Ttangent, Tstrain, eplp, 0.0, Eplp, eunn, funn, Enewn, 10);

					Trule = 10;

					return;
				}

				else							// on the compression envelope
				{
					envelope(x, Tstress, Ttangent, -1);

					TminStress = Tstress;
					TminStrain = Tstrain;

					return;
				}
			}

			else if (Trule == 14)
			{
				Tloading = 1;

				double er = Cstrain;		// Strain at reversal
				double fr = Cstress;		// stress at reversal
				double ea;					// strain where we rejoin rule 13
				double fa;					// stress where we rejoin rule 13
				double Enewn;			// reloading stiffness in tension;

				TReloadStrain = er;
				TReloadStress = fr;
				ea = T13Strain;
				fa = T13Stress;

				Enewn = funn/(eunn-epln);

				if (Tstrain > ea)
				{
					Ttangent = (fa - fr)/(ea - er);
					Tstress = Ttangent*(Tstrain - er) + fr;

					Trule = 15;

					return;
				}

				else if (Tstrain > eunn)		// Back on Rule 13
				{
					double Edum;

					Edum = (funn)/(eunn - T13Zero);
					
					calculate13Stress(Tstress, Ttangent, Tstrain, T13Zero, eunn, funn, Enewn);
					
					Trule = 13;
					return;
				}

				else							// Back on the compression envelope
				{
					envelope(x, Tstress, Ttangent, -1);

					TminStrain = Tstrain;
					TminStress = Tstress;

					return;
				}
			}

			else				// We are reversing post cracking
			{
				Tloading = 1;

				double er = Cstrain;		// Strain at reversal
				double fr = Cstress;		// stress at reversal
				double Enewn;			// reloading stiffness in tension;
				double Edum;
					
				if (eunn == 0 && funn == 0)
				{
					eunn = .05*epsc0;
					envelope(eunn/epsc0, funn, Edum, -1);
					TminStrain = eunn;
					TminStress = funn;
					Esecn = Ec*((funn/(Ec*epsc0)+0.57)/(eunn/epsc0+0.57));
					epln = eunn-funn/Esecn;
				}
					
				TReloadStrain = er;
				TReloadStress = fr;
				Enewn = funn/(eunn-epln);

				if (eunn/epsc0 >= xsp)
				{
					Tstress = 0.0;
					Ttangent = 0.0;
					Trule = 5;
					Tcracked = true;
					
					return;
				}


				if (Tstrain > eunn)
				{

					Edum = (funn) / (eunn - er);

					T13Zero = er;

					calculate13Stress(Tstress, Ttangent, Tstrain, er, eunn, funn, Enewn);

					Trule = 13;

					return;
				}

				else							// Back on the compression envelope
				{
					envelope(x, Tstress, Ttangent, -1);

					TminStrain = Tstrain;
					TminStress = Tstress;

					return;
				}
			}
		}
	}
}

double Concrete07::getStrain()
{
	return Tstrain;
}

double Concrete07::getStress()
{
	return Tstress;
}

double Concrete07::getTangent()
{
	return Ttangent;
}

int Concrete07::commitState()
{
	// History Variables
	CminStrain = TminStrain;
	CminStress = TminStress;
	CUnloadStrain = TUnloadStrain;
	CUnloadStress = TUnloadStress;
	CUnloadStiffness = TUnloadStiffness;
	CmaxStrain = TmaxStrain;
	CmaxStress = TmaxStress;
	CReloadStrain = TReloadStrain;
	CReloadStress = TReloadStress;
	CReloadStiffness = TReloadStiffness;
	C13Zero = T13Zero;
	C13Strain = T13Strain;
	C13Stress = T13Stress;
	Cloading = Tloading;
	Ccracked = Tcracked;
	Crule = Trule;

	// State Variables
	Cstrain = Tstrain;
	Cstress = Tstress;
	Ctangent = Ttangent;

	return 0;
}

int Concrete07::revertToLastCommit()
{
	// Reset History Variables to last committed state
	TminStrain = CminStrain;
	TminStress = CminStress;
	TUnloadStrain = CUnloadStrain;
	TUnloadStress = CUnloadStress;
	TUnloadStiffness = CUnloadStiffness;
	TmaxStrain = CmaxStrain;
	TmaxStress = CmaxStress;
	TReloadStrain = CReloadStrain;
	TReloadStress = CReloadStress;
	TReloadStiffness = CReloadStiffness;
	T13Zero = C13Zero;
	T13Strain = C13Strain;
	T13Stress = C13Stress;
	Tloading = Cloading;
	Tcracked = Ccracked;
	Trule = Crule;

	// Reset State Variables to last committed state
	Tstrain = Cstrain;
	Tstress = Cstress;
	Ttangent = Ctangent;

	return 0;
}

int Concrete07::revertToStart()
{

	// History Variables
	CminStrain = 0.0;
	CminStress = 0.0;
	CUnloadStrain = 0.0;
	CUnloadStress = 0.0;
	CUnloadStiffness = 0.0;
	CmaxStrain = 0.0;
	CmaxStress = 0.0;
	CReloadStrain = 0.0;
	CReloadStress = 0.0;
	CReloadStiffness = 0.0;
	C13Zero = 0.0;
	C13Strain = 0.0;
	C13Stress = 0.0;
	Cloading = 0.0;
	Ccracked = false;
	Crule = 0;

	TminStrain = 0.0;
	TminStress = 0.0;
	TUnloadStrain = 0.0;
	TUnloadStress = 0.0;
	TUnloadStiffness = 0.0;
	TmaxStrain = 0.0;
	TmaxStress = 0.0;
	TReloadStrain = 0.0;
	TReloadStress = 0.0;
	TReloadStiffness = 0.0;
	T13Zero = 0.0;
	T13Strain = 0.0;
	T13Stress = 0.0;
	Tloading = 0.0;
	Tcracked = false;
	Trule = 0;

	// State Variables
	Cstrain = 0.0;
	Cstress = 0.0;
	Ctangent = Ec;

	Tstrain = 0.0;
	Tstress = 0.0;
	Ttangent = Ec;

	return 0;
}

UniaxialMaterial* Concrete07::getCopy()
{
	Concrete07* theCopy = new Concrete07(this->getTag(), fpc, epsc0, Ec, fpt, epst0, xcrp, xcrn, r);

	// Converged history variables
	theCopy->CminStrain = CminStrain;
	theCopy->CminStress = CminStress;
	theCopy->CUnloadStrain = CUnloadStrain;
	theCopy->CUnloadStress = CUnloadStress;
	theCopy->CUnloadStiffness = CUnloadStiffness;
	theCopy->CmaxStrain = CmaxStrain;
	theCopy->CmaxStress = CmaxStress;
	theCopy->CReloadStrain = CReloadStrain;
	theCopy->CReloadStress = CReloadStress;
	theCopy->CReloadStiffness = CReloadStiffness;
	theCopy->C13Zero = C13Zero;
	theCopy->C13Strain = C13Strain;
	theCopy->C13Stress = C13Stress;
	theCopy->Cloading = Cloading;
	theCopy->Ccracked = Ccracked;
	theCopy->Crule = Crule;

	// Converged state variables
	theCopy->Cstress = Cstress;
	theCopy->Cstrain = Cstrain;
	theCopy->Ctangent = Ctangent;

	// Trial history variables
	theCopy->TminStrain = TminStrain;
	theCopy->TminStress = TminStress;
	theCopy->TUnloadStrain = TUnloadStrain;
	theCopy->TUnloadStress = TUnloadStress;
	theCopy->TUnloadStiffness = TUnloadStiffness;
	theCopy->TmaxStrain = TmaxStrain;
	theCopy->TmaxStress = TmaxStress;
	theCopy->TReloadStrain = TReloadStrain;
	theCopy->TReloadStress = TReloadStress;
	theCopy->TReloadStiffness = TReloadStiffness;
	theCopy->T13Zero = T13Zero;
	theCopy->Tloading = Tloading;
	theCopy->Tcracked = Tcracked;
	theCopy->Trule = Trule;

	// Trial state variables
	theCopy->Tstress = Tstress;
	theCopy->Tstrain = Tstrain;
	theCopy->Ttangent = Ttangent;

	return theCopy;
}

int Concrete07::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;
	static Vector data(31);
	data(0) = this->getTag();

	// Material Properties
	data(1) = fpc;
	data(2) = epsc0;
	data(3) = Ec;
	data(4) = fpt;
	data(5) = epst0;

	// Model Variables
	data(6) = xcrn;
	data(7) = xsp;
	data(8) = xcrp;
	data(9) = xcrk;
	data(10) = nn;
	data(11) = np;
	data(12) = r;

	// History Variables
	data(13) = CminStrain;
	data(14) = CminStress;
	data(15) = CUnloadStrain;
	data(16) = CUnloadStress;
	data(17) = CUnloadStiffness;
	data(18) = CmaxStrain;
	data(19) = CmaxStress;
	data(20) = CReloadStrain;
	data(21) = CReloadStress;
	data(22) = CReloadStiffness;
	data(23) = Cloading;
	if (Ccracked)
		data(24) = 1;
	else
		data(24) = 0;
	data(25) = Crule;

	// State Variables
	data(26) = Cstrain;
	data(27) = Cstress;
	data(28) = Ctangent;
	data(29) = C13Zero;
	data(30) = C13Strain;
	data(31) = C13Stress;

	// Data is only sent after convergence, so no trial variables
	// need to be sent through the data vector

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0)
		opserr << "Concrete07::sendSelf() - failed to send data\n";

	return res;
}

int Concrete07::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res = 0;
	static Vector data(28);
	res = theChannel.recvVector(this->getDbTag(), commitTag, data);

	if ( res < 0)
	{
		opserr << "Concrete07::recvSelf() - failed to receive data\n";
		this->setTag(0);
	}

	else
	{
		this->setTag(int(data(0)));

		// Material properties
		fpc = data(1);
		epsc0 = data(2);
		Ec = data(3);
		fpt = data(4);
		epst0 = data(5);

		// Model Variables
		xcrn = data(6);
		xsp = data(7);
		xcrp = data(8);
		xcrk = data(9);
		nn = data(10);
		np = data(11);
		r = data(12);

		// History Variables
		CminStrain = data(13);
		CminStress = data(14);
		CUnloadStrain = data(15);
		CUnloadStress = data(16);
		CUnloadStiffness = data(17);
		CmaxStrain = data(18);
		CmaxStress = data(19);
		CReloadStrain = data(20);
		CReloadStress = data(21);
		CReloadStiffness = data(22);
		Cloading = data(23);
		if (data(24) == 1)
			Ccracked = true;
		else
			Ccracked = false;
		Crule = data(25);

		// State Variables
		Cstrain = data(26);
		Cstress = data(27);
		Ctangent = data(28);
		C13Zero = data(29);
		C13Strain = data(30);
		C13Stress = data(31);

		Tstrain = Cstrain;
		Tstress = Cstress;
		Ttangent = Ctangent;
	}

	return res;
}

void Concrete07::Print(OPS_Stream &s, int flag)
{
	s << "Concrete07, tag: " << this->getTag() << endln;
	s << "  fpc: " << fpc << endln;
	s << "  epsc0: " << epsc0 << endln;
	s << "  fpt: " << fpt << endln;
	s << "  epst0: " << epst0 << endln;
	s << "  xsp: " << xsp << endln;
	s << "  xcrk: " << xcrk << endln;

	return;
}
