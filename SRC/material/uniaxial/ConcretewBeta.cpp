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

// Written: 2013, Yuan Lu & Marios Panagiotou (UC Berkeley) 

//consider making abstract class for materials that are sensitive to normal strain.

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <ConcretewBeta.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>
#include <float.h>

#include <string.h>
#include <stdlib.h>

#include <Information.h>
#include <MaterialResponse.h>

void *
	OPS_ConcretewBeta(void) {

		int total = OPS_GetNumRemainingInputArgs();

		if (total < 12) {
			opserr << "WARNING incorrect number of arguments\n";
			opserr << "Want: uniaxialMaterial ConcretewBeta $tag $fpc $ec0 $fcint $ecint $fcres $ecres $ft $ftint $etint $ftres $etres <-lambda $lambda> <-alpha $alpha> <-beta $bint $ebint $bres $ebres> <-E $E> <-conf $fcc ecc>\n";
			return 0;
		}

		int tag;
		double dData[11];
		double lambda = 0.5;
		double alpha = 1;
		double bData[4];
		bData[0] = 1;
		bData[1] = 0;
		bData[2] = 1;
		bData[3] = 0;
		double M = 0;
		double E0 = 0;
		double fcc = 0;
		double ecc = 0;

		int numData = 1;
		if (OPS_GetIntInput(&numData, &tag) != 0) {
			opserr << "WARNING invalid uniaxialMaterial Steel01 tag" << endln;
			return 0;
		}  

		numData = 11;
		if (OPS_GetDoubleInput(&numData, dData) != 0) {
			opserr << "WARNING invalid uniaxialMaterial Steel01 tag" << endln;
			return 0;
		}  

		total -= 12;

		//		char optionFlag[12];
		const char *optionFlag;
		while (total > 0) {
		  //			OPS_GetString(optionFlag,12);
		  optionFlag = OPS_GetString();
		  if (strcmp(optionFlag,"-beta") == 0) {
		    numData = 4;
		    if (OPS_GetDoubleInput(&numData, bData) != 0) {
		      opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument of -beta for tag " << tag << endln;
		      return 0;
		    }
		    total -= 5;
		  } else if (strcmp(optionFlag,"-lambda") == 0) {
				numData = 1;
				if (OPS_GetDoubleInput(&numData, &lambda) != 0) {
					opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument of -lambda for tag " << tag << endln;
					return 0;
				}
				total -= 2;
			} else if (strcmp(optionFlag,"-alpha") == 0) {
				numData = 1;
				if (OPS_GetDoubleInput(&numData, &alpha) != 0) {
					opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument of -alpha for tag " << tag << endln;
					return 0;
				}
				total -= 2;
			} else if (strcmp(optionFlag,"-M") == 0) {
				numData = 1;
				if (OPS_GetDoubleInput(&numData, &M) != 0) {
					opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument of -M for tag " << tag << endln;
					return 0;
				}
				total -= 2;
			} else if (strcmp(optionFlag,"-E") == 0) {
				numData = 1;
				if (OPS_GetDoubleInput(&numData, &E0) != 0) {
					opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument of -E for tag " << tag << endln;
					return 0;
				}
				total -= 2;
			} else if (strcmp(optionFlag,"-conf") == 0) {
				numData = 1;
				if (OPS_GetDoubleInput(&numData, &fcc) != 0) {
					opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument 1 of -conf for tag " << tag << endln;
					return 0;
				}
				if (OPS_GetDoubleInput(&numData, &ecc) != 0) {
					opserr << "WARNING invalid uniaxialMaterial ConcretewBeta argument 2 of -conf for tag " << tag << endln;
					return 0;
				}
				total -= 3;
			} else {
				opserr << "WARNING invalid uniaxialMaterial ConcretewBeta flag " << tag << endln;
				return 0;
			}
		}


		UniaxialMaterial *theMaterial = new ConcretewBeta(tag, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], lambda, alpha, bData[0], bData[1], bData[2], bData[3], M, E0, fcc, ecc);

		return theMaterial;

}

ConcretewBeta::ConcretewBeta(int tag, double fpc1, double ec01, double fcint1, double ecint1, double fcres1, 
	double ecres1, double fct1, double ftint1, double etint1, double ftres1, double etres1, 
	double lambda1, double alpha1, double bint1, double etbint1, double bres1, 
	double etbres1, double M1, double E01, double fcc1, double ecc1) 
	: UniaxialMaterial(tag, MAT_TAG_ConcretewBeta), 
	fpc(fpc1), ec0(ec01), fcint(fcint1), ecint(ecint1), fcres(fcres1), ecres(ecres1), 
	fct(fct1), ftint(ftint1), etint(etint1), ftres(ftres1), etres(etres1), 
	lambda(lambda1),  alpha(alpha1), bint(bint1), etbint(etbint1),
	bres(bres1), etbres(etbres1), M(M1), E0(E01), fcc(fcc1), ecc(ecc1)
{	
	if (fpc>0 || ec0>0 || fcint>0 || ecint>0 || fcres>0 || ecres>0) {
		opserr << "WARNING uniaxialMaterial ConcretewBeta tag " << tag << " has positive compression stress-strain values";
	} 
	if (ecint < ecres) {
		opserr << "WARNING uniaxialMaterial ConcretewBeta tag " << tag << " has ecint greater magnitude than ecres";		
	}
	if (fct1<0 || ftint<0 || etint<0 || ftres<0 || etres<0) {
		opserr << "WARNING uniaxialMaterial ConcretewBeta tag " << tag << " has negative tension stress-strain values";
	}
	if (lambda<0 || lambda>1) {
		opserr << "WARNING uniaxialMaterial ConcretewBeta tag " << tag << " has invalid lambda value";
	}
	if (alpha<0) {
		opserr << "WARNING uniaxialMaterial ConcretewBeta tag " << tag << " has invalid alpha value";
	}
	if (bint<0 || etbint<0 || bres<0 || etbres<0) {
		opserr << "WARNING uniaxialMaterial ConcretewBeta tag " << tag << " has negative beta values";
	}
	if (fcc == 0) {
		fcc = fpc;
		ecc = ec0;
	} else if (fcc > fpc) {
		opserr << "WARNING uniaxialMaterial ConcretewBeta tag " << tag << " has fcc smaller than fpc: material will ignoring confinement";
		fcc = fpc;
		ecc = ec0;
	} else if (ecc > ec0) {
		opserr << "WARNING uniaxialMaterial ConcretewBeta tag " << tag << " has ecc smaller than ec0: material will ignoring confinement";
		fcc = fpc;
		ecc = ec0;
	}
	if (M < 0) {
		opserr << "WARNING uniaxialMaterial ConcretewBeta tag " << tag << " has invalid M value";
	}
	if (E0 == 0) {
		E0 = 2*(fpc/ec0);
	} else if (E0 < (fpc/ec0)) {
		opserr << "WARNING uniaxialMaterial ConcretewBeta tag " << tag << " too small E0 value - setting to (fpc/ec0)";
		E0 = (fpc/ec0);
	} else if (E0 > (2*fpc/ec0)) {
		opserr << "WARNING uniaxialMaterial ConcretewBeta tag " << tag << " too large E0 value - setting to (2*fpc/ec0)";
		E0 = 2*(fpc/ec0);
	}

	this->updateStoredValues();
	this->revertToStart();
}	

ConcretewBeta::ConcretewBeta()
	: UniaxialMaterial(0, MAT_TAG_ConcretewBeta), 
	fpc(0), ec0(0), E0(0), fcint(0), ecint(0), fcres(0), 
	ecres(0), fct(0), ftint(0), etint(0), ftres(0), etres(0), 
	lambda(0),  alpha(0), bint(0), etbint(0),
	bres(0), etbres(0), M(0)
{
	// You probably shouldn't be calling this anyways. Just seems like a bad idea.
	this->updateStoredValues();
	this->revertToStart();
}

ConcretewBeta::~ConcretewBeta() 
{
	// nothing for now
}

int ConcretewBeta::setTrialStrain(double strain, double strainRate) {
	return setTrialStrainwBeta(strain, 0.0, strainRate);
}

// extra thing, only call if you know what you're doing.
int ConcretewBeta::setTrialStrainwBeta(double newStrain, double et, double strainRate) {
	//rever to last value.
	this->revertToLastCommit();

	double newStress, newStressPure, newTangent;
	// need to set beta.
	double beta = computeBeta(newStrain, et);
	//double beta = 0.7;
	if ( fabs(newStrain - Tstrain) < DBL_EPSILON) {
		return 0;
	}
	setValues(newStrain, beta, newStress, newStressPure, newTangent);
	Ttangent = newTangent;
	Tstress = newStress;
	Tstrain = newStrain;
	Tbeta = beta;
	if (Tstrain >= CMaxStrainTens) {
		TMaxStrainTens = Tstrain;
		TMaxStressTens = Tstress;
	} else if (Tstrain <= CMaxStrainCompr) {
		TMaxStrainCompr = Tstrain;
		TMaxStressCompr = Tstress;
		TMaxStressComprPure = newStressPure;
	}
	return 0;
}

double ConcretewBeta::getStrain(void)
{
	return Tstrain;
}

double ConcretewBeta::getStress(void) 
{
	return Tstress;
}

double ConcretewBeta::getTangent(void)
{
	return Ttangent;
}

double ConcretewBeta::getBeta(void)
{
	return Tbeta;
}

double ConcretewBeta::getInitialTangent(void) 
{
	return E0;
}

int ConcretewBeta::commitState(void) 
{	
	CMaxStrainCompr = TMaxStrainCompr;
	CMaxStressCompr = TMaxStressCompr;
	CMaxStressComprPure = TMaxStressComprPure;
	CMaxStrainTens = TMaxStrainTens;
	CMaxStressTens = TMaxStressTens;

	Cstrain = Tstrain;
	Cstress = Tstress;
	Ctangent = Ttangent;
	Cbeta = Tbeta;

	return 0;
}

int ConcretewBeta::revertToLastCommit(void)
{
	TMaxStrainCompr = CMaxStrainCompr;
	TMaxStressCompr = CMaxStressCompr;
	TMaxStressComprPure = CMaxStressComprPure;
	TMaxStrainTens = CMaxStrainTens;
	TMaxStressTens = CMaxStressTens;

	Tstrain = Cstrain;
	Tstress = Cstress;
	Ttangent = Ctangent;
	Tbeta = Cbeta;

	return 0;
}

int ConcretewBeta::revertToStart(void) 
{
	CMaxStrainCompr = 0.0;
	CMaxStressCompr = 0.0;
	CMaxStressComprPure = 0.0;
	CMaxStrainTens = 0.0;
	CMaxStressTens = 0.0;

	Cstrain = 0.0;
	Cstress = 0.0;
	Ctangent = E0;
	Cbeta = 1;

	return this->revertToLastCommit();
}

UniaxialMaterial * ConcretewBeta::getCopy(void) 
{
	ConcretewBeta* theCopy = new ConcretewBeta(this->getTag(), this->fpc, this->ec0, this->fcint, this->ecint, 
		this->fcres, this->ecres, this->fct, this->ftint, this->etint, this->ftres, this->etres, this->lambda, 
		this->alpha, this->bint, this->etbint, this->bres, this->etbres, this->M, this->E0, this->fcc, this->ecc);


	// Converged history variables
	theCopy->CMaxStrainCompr = CMaxStrainCompr;
	theCopy->CMaxStressCompr = CMaxStressCompr;
	theCopy->CMaxStressComprPure = CMaxStressComprPure;
	theCopy->CMaxStrainTens = CMaxStrainTens;
	theCopy->CMaxStressTens = CMaxStressTens;

	// Converged state variables
	theCopy->Cstrain = Cstrain;
	theCopy->Cstress = Cstress;
	theCopy->Ctangent = Ctangent;
	theCopy->Cbeta = Cbeta;

	theCopy->revertToLastCommit();

	return theCopy;
}

int ConcretewBeta::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;
	static Vector data(40);
	data(0) = this->getTag();

	data(1) = fpc;
	data(2) = ec0;
	data(3) = fcint;
	data(4) = ecint;
	data(5) = fcres;
	data(6) = ecres;
	data(7) = fct;
	data(8) = ftres;
	data(9) = etres;
	data(10) = ftres;
	data(11) = etres;
	data(12) = E0;
	data(13) = fcc;
	data(14) = ecc;

	// History variables from last converged state
	data(20) = CMaxStrainCompr;
	data(21) = CMaxStressCompr;
	data(22) = CMaxStressComprPure;
	data(23) = CMaxStrainTens;
	data(24) = CMaxStressTens;

	// State variables from last converged state
	data(25) = Cstrain;
	data(26) = Cstress;
	data(27) = Ctangent;
	data(28) = Cbeta;

	// other values:
	data(29) = lambda;
	data(30) = alpha;
	data(31) = bint;
	data(32) = etbint;
	data(33) = bres;
	data(34) = etbres;
	data(35) = M;

	// Data is only sent after convergence, so no trial variables
	// need to be sent through data vector

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) 
		opserr << "ConcretewBeta::sendSelf() - failed to send data\n";

	return res;
}

int ConcretewBeta::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) 
{
	int res = 0;
	static Vector data(40);
	res = theChannel.recvVector(this->getDbTag(), commitTag, data);

	if (res < 0) {
		opserr << "ConcretewBeta::recvSelf() - failed to receive data\n";
		this->setTag(0);      
	} else {
		this->setTag(int(data(0)));

		// Material properties 
		fpc = data(1);
		ec0 = data(2);
		fcint = data(3);
		ecint = data(4);
		fcres = data(5);
		ecres = data(6);
		fct = data(7);
		ftint = data(8);
		etint = data(9);
		ftres = data(10);
		etres = data(11);
		E0 = data(12);
		fcc = data(13);
		ecc = data(14);

		// Extra saved stuff
		this->updateStoredValues();

		// History variables from last converged state
		CMaxStrainCompr = data(20);
		CMaxStressCompr = data(21);
		CMaxStressComprPure = data(22);
		CMaxStrainTens = data(23);
		CMaxStressTens = data(24);

		// State variables from last converged state
		Cstrain = data(25);
		Cstress = data(26);
		Ctangent = data(27);
		Cbeta = data(28);

		// other values
		lambda = data(29);
		alpha = data(30);
		bint = data(31);
		etbint = data(32);
		bres = data(33);
		etbres = data(34);
		M = data(35);

		// Set trial state variables
		this->revertToLastCommit();
	}

	return res;
}

void ConcretewBeta::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"ConcretewBeta\", ";
		s << "\"fpc\": " << fpc << ", ";
		s << "\"ec0\": " << ec0 << ", ";
		s << "\"fcint\": " << fcint << ", ";
		s << "\"ecint\": " << ecint << ", ";
		s << "\"fcres\": " << fcres << ", ";
		s << "\"ecres\": " << ecres << ", ";
		s << "\"ft\": " << fct << ", ";
		s << "\"ftint\": " << ftint << ", ";
		s << "\"etint\": " << etint << ", ";
		s << "\"ftres\": " << ftres << ", ";
		s << "\"etres\": " << etres << ", ";
		s << "\"lambda\": " << lambda << ", ";
		s << "\"alpha\": " << alpha << ", ";
		s << "\"bint\": " << bint << ", ";
		s << "\"ebint\": " << etbint << ", ";
		s << "\"bres\": " << bres << ", ";
		s << "\"ebres\": " << etbres << ", ";
		s << "\"E\": " << E0 << ", ";
		s << "\"fcc\": " << fcc << ", ";
		s << "\"ecc\": " << ecc << "}";
	} else {
		s << "ConcretewBeta, tag: " << this->getTag() << endln;
	}
}


int ConcretewBeta::setValues(double newStrain, double beta, double & newStress, double & newStressPure, double & newTangent)
{	
	if (E0 == 0) {
		// case where material is just 0.
		newStressPure = 0;
		newStress = 0;
		newTangent = 0;
		return 0;
	}

	if (Tstress >= 0 && beta != 1.0) {
		// no beta in tension range.
		beta = 1.0;
	}

	double de = newStrain - Tstrain;
	if (de > 0) {
		// going in tension direction
		if (newStrain > CMaxStrainTens) {
			// loading in tension
			if (newStrain <= et0) {
				newStressPure = E0*newStrain;
				newStress = newStressPure;
				newTangent = E0;
				//} else {
				//	newStressPure = (fct - ftres)*pow(et0/newStrain,bexp) + ftres;
				//	newStress = newStressPure;
				//	newTangent = (newStress - Tstress)/(newStrain - Tstrain);
				//}
			} else if (M == 0) {
				if (newStrain <= etint) {
					newTangent = (fct - ftint)/(et0 - etint);
					newStressPure = newTangent*(newStrain - etint) + ftint;
					newStress = newStressPure;
				} else if (newStrain <= etres) {
					newTangent = (ftint - ftres)/(etint - etres);
					newStressPure = newTangent*(newStrain - etres) + ftres;
					newStress = newStressPure;
				} else {
					newStressPure = ftres;
					newStress = newStressPure;
					newTangent = 0;
				}
			} else {
				//newStressPure = fct/(1+sqrt(3.6*M*(newStrain-et0))); // Bentz tension stiffening.
				if (lambdaM == 1.0) {
					newStressPure = fct;
					newStress = newStressPure;
					newTangent = 0;
				} else {
					newStressPure = fct*((1-M)*exp(-lambdaM*(newStrain-et0))+M); // Stevens-Collins tension stiffening
					newStress = newStressPure;
					newTangent = -lambdaM*fct*((1-M)*exp(-lambdaM*(newStrain-et0)));
				}
			}
			//} else if ((Tstress < 0) && (eaftc!=0) && (CMaxStrainTens>(fct/E0)) && (alphafct/eaftc)*Tstrain >= Tstress) {
			//	// reloading in tension while in "unloading in tension" range...
			//	newStress = Tstrain + E0*de;
			//	newTangent = E0;
			//	if (newStress > 0) {
			//		double crossPoint = Tstrain + (0-Tstress)/E0;
			//		double tang2 = (CMaxStressTens)/(CMaxStrainTens - crossPoint);
			//		newStress = Tstress + E0*(crossPoint-Tstrain) + tang2*(newStrain-crossPoint);
			//		newTangent = (newStress-Tstress)/(newStrain-Tstrain);
			//	}
			//	newStressPure = newStress;
		} else if (Tstress < 0) {
			// unloading in compression
			double E1;
			if (Tstrain != 0) {
				E1 = (Tstress/Tstrain);
			} else {
				E1 = 0;
			}
			newTangent = ((1 - lambda)*E0 + lambda*E1);
			newStressPure = newTangent*de + Tstress;
			newStress = newStressPure;
			if (E1 > E0 || E1 < 0) {
				// go for max tension directly
				newTangent = (CMaxStressTens - Tstress)/(CMaxStrainTens - Tstrain);
				newStressPure = newTangent*de + Tstress;
				newStress = newStressPure;
			} else if (newStress > 0) {
				double zeroStrain = newStrain - newStress/newTangent;
				newTangent = (CMaxStressTens)/(CMaxStrainTens - zeroStrain);
				newStressPure = newTangent*(newStrain - zeroStrain);
				newStress = newStressPure;
			}
		} else {
			// reloading in tension
			// To handle random numerical jumps due to rounding errors.
			newTangent = (CMaxStressTens - Tstress)/(CMaxStrainTens - Tstrain);
			double tangOption = (Tstress)/(Tstrain);
			if (newStrain < 0 && tangOption > newTangent) {
				// force it to cross (0, 0)
				newTangent = tangOption;
			} 
			//if (newTangent < 0) {
			//	// Numerical problem due to rounding. 
			//	newTangent = 0;
			//	newStressPure = Tstress;
			//	newStress = Tstress;
			//} else {
			newStressPure = newTangent*de + Tstress;
			newStress = newStressPure;
			//}
		}
	} else if (de < 0) {
		// going in compression direction
		if ((newStrain - CMaxStrainCompr) <= DBL_EPSILON) {
			// loading in Compression
			if (newStrain >= eaftc && CMaxStrainTens > 0) {
				// reloading due to tension
				newTangent = (alphafct - Tstress)/(eaftc - Tstrain);
				newStressPure = newTangent*de + Tstress;
				newTangent = (beta*alphafct - Tstress)/(eaftc - Tstrain);
				newStress = newTangent*de + Tstress;
			} else if (newStrain >= ec0) {
				// loading normally by equation
				double temp2 = (fpc/(ec0*ec0)-E0/ec0);
				newStressPure = E0*newStrain + temp2*(newStrain*newStrain);
				// newStressPure = fpc*((2*newStrain/ec0) - pow(newStrain/ec0,2));
				// newTangent = beta*newTangentPure;
				newStress = beta*newStressPure;
				newTangent = beta*(E0 + 2*temp2*newStrain);
			} else if (newStrain >= ecc) {
				// fpc -> fcc branch for confinement
				double temp2 = (fpc-fcc)/pow(ec0-ecc,3)*pow(newStrain-ecc,2);
				newStressPure = temp2*(newStrain-ecc)+fcc;
				//newStressPure = aLoad*pow(newStrain,3)+bLoad*pow(newStrain,2)+cLoad*newStrain+dLoad;
				newStress = beta*newStressPure;
				newTangent = 3*temp2;
			} else if (newStrain >= ecint) {
				// linear softening
				newStressPure = ElinearsoftcP1*(newStrain - ecint) + fcint;
				newStress = beta*newStressPure;
				newTangent = beta*ElinearsoftcP1;
			} else if (newStrain >= ecres) {
				// linear softening
				newStressPure = ElinearsoftcP2*(newStrain - ecres) + fcres;
				newStress = beta*newStressPure;
				newTangent = beta*ElinearsoftcP2;
			} else {
				newStressPure = fcres;
				newStress = beta*newStressPure;
				//if ( fabs(newStrain - Tstrain) < DBL_EPSILON) {
				//	newTangent = 0;
				//	newStressPure = Tstress/Tbeta;
				//	newStress = beta*newStressPure;
				//} else {
				newTangent = (newStress - Tstress)/(newStrain - Tstrain);
				//}
			}
		} else if ((E0*de+Tstress) >= 0) {
			// unload in tension to stress 0
			newTangent = E0;
			newStress = E0*de+Tstress;
			newStressPure = newStress;
		} else {
			// reloading in Compression
			double unloadStress = 0;
			double unloadStrain = 0;
			if (Tstress > 0) {
				unloadStress = Tstress;
				unloadStrain = Tstress/E0;
				Tstrain = Tstrain - unloadStrain;
				if (Tstrain < newStrain) {
					Tstrain = newStrain;
					unloadStrain = (Tstrain - newStrain);
					unloadStress = E0*unloadStrain;
				}
				Tstress = Tstress - unloadStress;
			}
			if (Tstrain >= eaftc) {
				// aim for either that one point or max compression
				double E1a = (alphafct - Tstress)/(eaftc - Tstrain);
				double E1 = (beta*alphafct - Tstress)/(eaftc - Tstrain);
				double E2, E2a;
				double maxaim;
				if (beta*CMaxStressComprPure > CMaxStressCompr) {
					E2a = (CMaxStressComprPure - Tstress)/(CMaxStrainCompr - Tstrain);
					maxaim = beta*CMaxStressComprPure;
					E2 = (maxaim - Tstress)/(CMaxStrainCompr - Tstrain);
				} else {
					maxaim = CMaxStressCompr;
					E2 = (maxaim - Tstress)/(CMaxStrainCompr - Tstrain);
					E2a = (maxaim/beta - Tstress)/(CMaxStrainCompr - Tstrain);
				}

				//E2a = (CMaxStressComprPure - Tstress)/(CMaxStrainCompr - Tstrain);
				//maxaim = beta*CMaxStressComprPure;
				//E2 = (maxaim - Tstress)/(CMaxStrainCompr - Tstrain);

				//maxaim = CMaxStressCompr;
				//double E2 = (maxaim - Tstress)/(CMaxStrainCompr - Tstrain);
				//if (newStrain > 0 && Tstrain < 0) {
				//	if (E1 > E2 && CMaxStrainCompr <= ec0) {
				//		newStressPure = E2a*de + Tstress;
				//		newStress = E2*(newStrain) + E2a*(-Tstrain) + Tstress;
				//		newTangent = (newStress - Tstress)/(newStrain - Tstrain);
				//	} else {
				//		newStressPure = E1*de + Tstress;
				//		newStress = E1*(newStrain) + E1a*(-Tstrain) + Tstress;
				//	}
				//} else {
				if (E1a > E2a && CMaxStrainCompr <= ec0) {
					newTangent = (maxaim - (Tstress+unloadStress))/(CMaxStrainCompr - (Tstrain+unloadStrain));
					newStressPure = E2a*(newStrain - Tstrain) + Tstress;
					newStress = E2*(newStrain - Tstrain) + Tstress;
				} else {
					newTangent = (beta*alphafct - (Tstress+unloadStress))/(eaftc - (Tstrain+unloadStrain));
					newStressPure = E1a*(newStrain - Tstrain) + Tstress;
					newStress = E1*(newStrain - Tstrain) + Tstress;
				}
				//}
			} else {
				// aim for max compression
				double maxaimA, maxaim, E2, E2a;
				if (beta*CMaxStressComprPure > CMaxStressCompr) {
					maxaim = beta*CMaxStressComprPure;
					maxaimA = CMaxStressComprPure;
				} else {
					maxaim = CMaxStressCompr;
					maxaimA = CMaxStressCompr/beta;
				}
				maxaim = beta*CMaxStressComprPure;
				maxaimA = CMaxStressComprPure;

				//Tstress = Tstress + unloadStress;
				//Tstrain = Tstrain + unloadStrain;
				//if (Tstrain > eaftc) {
				//	newStress = (maxaim - beta*alphafct)/(CMaxStrainCompr - eaftc)*(newStrain - eaftc) + beta*alphafct;
				//	newStressPure = (maxaimA - alphafct)/(CMaxStrainCompr - eaftc)*(newStrain - eaftc) + alphafct;
				//	newTangent = (newStress - Tstress)/(newStrain - Tstrain);
				//} else {
				E2 = (maxaim - (Tstress + unloadStress))/(CMaxStrainCompr - (Tstrain + unloadStrain));
				E2a = (maxaimA - (Tstress + unloadStress))/(CMaxStrainCompr - (Tstrain + unloadStrain));
				newTangent = E2;
				newStressPure = E2a*(newStrain - Tstrain) + Tstress;
				newStress = E2*(newStrain - Tstrain) + Tstress;
				//}
			}
		} 
	} else {
		if (Tstress <= 0) {
			// compression.
			newStressPure = Tstress/Tbeta;
			newStress = Tstress;
			newTangent = Ttangent;
		} else {
			// tension
			newStressPure = Tstress;
			newStress = Tstress;
			newTangent = Ttangent;
		}
	}

	//if (newTangent > E0) {
	//	int x = 1;
	//} else if (newTangent < -E0) {
	//	int x = 1;
	//}
	return 0;
}

double ConcretewBeta::computeBeta(double newStrain, double et) {
	if (et <= 0.0) {
		return 1.0;
	} else if (et <= etbint) {
		return (1.0 - bint)/(-etbint)*(et - etbint) + bint;;
	} else if (et <= etbres) {
		return (bint - bres)/(etbint - etbres)*(et - etbres) + bres;
	} else {
		return bres;
	}

	//if (et <= 0.0) {
	//	return 1.0;
	//} 
	//double beta = 1/(0.8-0.34*(et/ec0));
	//if (beta <= 1.) {
	//	return beta;
	//} else {
	//	return 1.0;
	//}
}

void ConcretewBeta::updateStoredValues() 
{
	if (ec0 == 0) {	
		et0 = 0;
		alphafct = 0;
		eaftc = 0;
		ElinearsoftcP1 = 0;
		ElinearsoftcP2 = 0;
		lambdaM = 1.0;
		//aLoad = 0;
		//bLoad = 0;
		//cLoad = 0;
		//dLoad = fpc;
	} else {
		et0 = fct/E0;
		alphafct = -alpha*fct;
		double e1 = (-E0 + sqrt(E0*E0 - 4*(fpc-E0*ec0)/(ec0*ec0)*(-alphafct)))/(fpc-E0*ec0)*0.5*(ec0*ec0);
		double e2 = (-E0 - sqrt(E0*E0 - 4*(fpc-E0*ec0)/(ec0*ec0)*(-alphafct)))/(fpc-E0*ec0)*0.5*(ec0*ec0);
		if (e1 > e2 && e1 <= 0) {
			eaftc = e1;
		} else {
			eaftc = e2;
		}
		if (ec0 != ecint) {
			ElinearsoftcP1 = (fcc - fcint)/(ecc - ecint);
		} else {
			ElinearsoftcP1 = 0;
		}
		if (ecint != ecres) {
			ElinearsoftcP2 = (fcint - fcres)/(ecint - ecres);
		} else {
			ElinearsoftcP2 = 0;
		}
		if (M > 0.) {
			lambdaM = 540./sqrt(M);
			//if (lambdaM > 1000.) {
			//	lambdaM = 1000.;
			//}
		} else {
			lambdaM = 1.0;
		}
		//if ((fcc == fpc) || (ec0 == ecc)) {
		//	aLoad = 0;
		//	bLoad = 0;
		//	cLoad = 0;
		//	dLoad = fpc;
		//} else {
		//	double temp = pow(ec0-ecc,3);
		//	aLoad = -2*(fpc-fcc)/temp;
		//	bLoad = 3*(ec0+ecc)*(fpc-fcc)/temp;
		//	cLoad = -(6*ec0*ecc*(fpc-fcc))/temp;
		//	dLoad = (fcc*pow(ec0,3)-3*fcc*ec0*ec0*ecc + 3*fpc*ec0*ecc*ecc - fpc*pow(ecc,3))/temp;
		//}
	}

}

Response* 
	ConcretewBeta::setResponse(const char **argv, int argc,
	OPS_Stream &theOutput)
{
	Response *theResponse = 0;

	theOutput.tag("UniaxialMaterialOutput");
	theOutput.attr("matType", this->getClassType());
	theOutput.attr("matTag", this->getTag());

	// stress
	if (strcmp(argv[0],"stress") == 0) {
		theOutput.tag("ResponseType", "sigma11");
		theResponse =  new MaterialResponse(this, 1, this->getStress());
	}  
	// tangent
	else if (strcmp(argv[0],"tangent") == 0) {
		theOutput.tag("ResponseType", "C11");
		theResponse =  new MaterialResponse(this, 2, this->getTangent());
	}

	// strain
	else if (strcmp(argv[0],"strain") == 0) {
		theOutput.tag("ResponseType", "eps11");
		theResponse =  new MaterialResponse(this, 3, this->getStrain());
	}

	// strain
	else if ((strcmp(argv[0],"stressStrain") == 0) || 
		(strcmp(argv[0],"stressANDstrain") == 0) ||
		(strcmp(argv[0],"stressAndStrain") == 0)) {
			theOutput.tag("ResponseType", "sig11");
			theOutput.tag("ResponseType", "eps11");
			theResponse =  new MaterialResponse(this, 4, Vector(2));
	}

	else if ((strcmp(argv[0],"stressStrainTangent") == 0) || 
		(strcmp(argv[0],"stressANDstrainANDtangent") == 0)) {
			theOutput.tag("ResponseType", "sig11");
			theOutput.tag("ResponseType", "eps11");
			theOutput.tag("ResponseType", "C11");
			theResponse =  new MaterialResponse(this, 5, Vector(3));
	}

	// beta
	else if (strstr(argv[0],"beta") != 0) {
		theOutput.tag("ResponseType", "beta");
		theResponse =  new MaterialResponse(this, 6, this->getBeta());
	}

	// stress sensitivity for local sensitivity recorder purpose.  Quan 2009
	// limit:  no more than 10000 random variables/sensitivity parameters
	else if (strstr(argv[0],"stressSensitivity") != 0) {
		char *token = strtok((char *) argv[0], " ");
		if (token != NULL) token = strtok(NULL, " ");
		int gradient = atoi(token);
		theOutput.tag("ResponseType", "sigsens11");
		theResponse =  new MaterialResponse(this, gradient+10000, this->getStress());
	}
	// strain sensivitiy
	else if (strstr(argv[0],"strainSensitivity") != 0) {
		char *token = strtok((char *) argv[0], " ");
		if (token != NULL) token = strtok(NULL, " ");
		int gradient = atoi(token);
		theOutput.tag("ResponseType", "epssens11");
		theResponse =  new MaterialResponse(this, gradient+20000, this->getStrain());
	}

	theOutput.endTag();
	return theResponse;

}

int 
	ConcretewBeta::getResponse(int responseID, Information &matInfo)
{
	static Vector stressStrain(2);
	static Vector stressStrainTangent(3);
	// each subclass must implement its own stuff   

	// added for sensitivity recorder. Quan 2009
	if ((responseID>10000)&&(responseID<20000)){
		matInfo.setDouble(this->getStressSensitivity(responseID-10000,false));
		return 0;
	}
	else if (responseID>20000){
		matInfo.setDouble(this->getStrainSensitivity(responseID-20000));
		return 0;
	}

	switch (responseID) {
	case 1:
		matInfo.setDouble(this->getStress());
		return 0;

	case 2:
		matInfo.setDouble(this->getTangent());
		return 0;      

	case 3:
		matInfo.setDouble(this->getStrain());
		return 0;      

	case 4:
		stressStrain(0) = this->getStress();
		stressStrain(1) = this->getStrain();
		matInfo.setVector(stressStrain);
		return 0;

	case 5:
		stressStrainTangent(0) = this->getStress();
		stressStrainTangent(1) = this->getStrain();
		stressStrainTangent(2) = this->getTangent();
		matInfo.setVector(stressStrainTangent);
		return 0;

	case 6:
		matInfo.setDouble(this->getBeta());
		return 0;  

	default:      
		return -1;
	}
}
