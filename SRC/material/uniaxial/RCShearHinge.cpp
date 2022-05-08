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

// $Revision: 1.0 $
// $Date: 2022/05/02 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/RCShearHinge.cpp,v $

// Written: Amir Reza Tabkhi Wayghan, MASc, Structural Engineering Graduate, Carleton University/rCon Engineers Inc.
//          & Vahid Sadeghian, PhD, Assistant Professor, Carleton University
// Created: May, 2022
// Revision: A
//
//
// Description: This file contains the implementation for the RCShearHinge model developed based on the following journal papers & thesis:
// 				Tabkhi, A.R. and Sadeghian, V. (accepted) “A Shear Hinge Model for Analysis of Reinforced Concrete Columns,” ACI Structural Journal.
//				Tabkhi, A.R. and Sadeghian, V. (2021) “A Shear Hinge Model for Analysis of Reinforced Concrete Beams,” ACI Structural Journal, Vol. 118, No. 6, pp. 279-291.
//				Tabkhi, A.R. (2021) "Development of Shear Plastic Hinge Models for Analysis of Reinforced Concrete Members," MASc Thesis, Carleton University, 2021.

#include <math.h>
#include <string.h>

#include <elementAPI.h>
#include <RCShearHinge.h>
#include <Vector.h>
#include <Channel.h>
#include <MaterialResponse.h>

#include <OPS_Globals.h>
#define PI 3.14159265359
static int numRCShearHingeMaterials = 0;

void*
OPS_RCShearHinge()
{
	if (numRCShearHingeMaterials == 0) {
		numRCShearHingeMaterials++;
		opserr << "RCShearHingeMaterial, Developed by Amir Reza Tabkhi & Vahid Sadeghian \n";
	}

	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	int  tag;
	double isDeep = 0;

	double dData[20];
	int numData = 1;
	// Check tag
	if (OPS_GetIntInput(&numData, &tag) != 0) {
		opserr << "WARNING invalid uniaxialMaterial  RCShearHinge tag" << endln;
		return 0;
	}

	numData = OPS_GetNumRemainingInputArgs();		

	if (numData == 20) {
		isDeep = 1;
	 }

	if (numData == 19) {
		isDeep = 0;
	}

	if (numData != 19 && numData != 20) {
		opserr << "Invalid Args want: uniaxialMaterial RCShearHinge tag b h d pRat fc Ec eps0 fyt fyl Es Ast Asl Asc s a alpha cover forcecf dispcf <isDeep?>\n ";
		return 0;
	}

	if (OPS_GetDoubleInput(&numData, dData) != 0) {            
		opserr << "Invalid Args want: uniaxialMaterial RCShearHinge tag b h d pRat fc Ec eps0 fyt fyl Es Ast Asl Asc s a alpha cover forcecf dispcf <isDeep?>\n ";
		return 0;
	}
	if (numData == 21)
	{
		if (OPS_GetIntInput(&numData, &tag) != 0 && isDeep != 0) {
			opserr << "WARNING invalid isTag- 0/1 integer is required" << endln;
			return 0;
		}
	}
	// Parsing was successful, allocate the material with zero index
	theMaterial = new RCShearHinge(tag,
		dData[0], dData[1], dData[2], dData[3],
		dData[4], dData[5], dData[6], dData[7],
		dData[8], dData[9], dData[10], dData[11],
		dData[12], dData[13], dData[14], dData[15],
		dData[16], dData[17], dData[18], isDeep);

	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type RCShearHinge \n";
		return 0;
	}

	return theMaterial;
}



RCShearHinge::RCShearHinge(int tag, double b, double h, double d, double pRat,
	double fc, double Ec, double eps0, double fyt, double fyl, double Es,
	double Ast, double Asl, double Asc, double s, double a, double alpha, double cover, double forcecf, double dispcf, double isDeep) :UniaxialMaterial(tag, MAT_TAG_RCShearHinge)
{
	// note that pRat is considered as positive for compression in this code. 
	// IsDeep is 0 (ignore=default) or 1 (consider)


	// general factors & parameters
	double dv = std::max(0.9 * d, 0.72 * h);
	double exay = Ast * fyt * (dv / s) / (pow(Ast * fyt / (0.1 * fc * b * s), 0.23));

	// ULTIMATE POINT
	double k1 = 750 * (1 + alpha) / (Asl * Es);
	double k1p = -375 * (pRat * b * h * fc) / (Asl * Es);
	double k2 = 0.4 * sqrt(fc) * b * dv;
	double k3 = (1.73 - 0.2 * k1p) * exay;
	double k4 = 1 + 0.2 * exay * k1;
	double k15 = Ast * fyt * dv / s;

	// Calculation of Vu
	double vu1 = (k1 * k3 - k4 * (1 + k1p)) / (2 * k1 * k4) + sqrt(pow((k1 * k3 - k4 * (1 + k1p)) / (2 * k1 * k4), 2) + (k2 + k3 * (1 + k1p)) / (k1 * k4));
	double vu2 = sqrt(pow(alpha * k15, 2) + 2 * k15 * Asl * fyl + k15 * pRat * b * h * fc) - alpha * k15 + 0.05 * b * dv * sqrt(fc);
	vu = std::min(vu1, vu2);

	// Calculation of exu
	double exu_ten = (k1 * vu + k1p) / 750;
	double k7p = Asc * exu_ten * Es;
	double c = std::min((1 + alpha) * vu - 0.5 * pRat * b * h * fc, Asl * fyl) + pRat * b * h * fc;
	double Xu = (sqrt(k7p * k7p * eps0 * eps0 / 4 + k7p * eps0 * (b * (h - d) * fc * exu_ten + c * eps0 / 2) + c * eps0 * (b * d * fc * exu_ten + c * eps0 / 4)) - (k7p + c) * eps0 / 2) / (b * fc * exu_ten);
	double k5 = (h / 2 - Xu) / (d - Xu);
	double exu = std::max((k1 * vu + k1p) * k5 / 750, -1 * pRat * fc / (0.5 * Ec + Asl * Es / (b * d)));

	// Calculation of thetau
	double pzfyfc = Ast * fyt / (fc * b * s);
	double k6 = pzfyfc < 0.1 ? pow(10 * pzfyfc, 0.2) : 1;

	double thetau = (29 + 7000 * exu) * k6;

	// Calculation of e2u
	double thetau_ten = thetau / (1 - 2 * Xu / (3 * h));
	double thetau_comp = thetau_ten / 3;
	double fc2u_ten = vu / (b * dv) * (tan(thetau_ten * PI / 180.0) + 1 / (tan(thetau_ten * PI / 180.0)));
	double fc2u_comp = vu / (b * dv) * (tan(thetau_comp * PI / 180.0) + 1 / (tan(thetau_comp * PI / 180.0)));
	double fc2u = (fc2u_comp * Xu + fc2u_ten * (h - Xu)) / h;
	double delta = sqrt(1 + pow(tan(2 * thetau * PI / 180.0), 2));
	double miu = 85 * (delta + 1) / (delta - 1) * fc2u / fc;
	double k8 = miu * eps0 - 1 + sqrt(pow((miu * eps0 - 1), 2) + 4 * miu * exu * delta / (delta + 1) + 0.8 * fc2u / fc);
	double e2u = std::max(k8, 1 - sqrt(1 - fc2u / fc)) * eps0;

	// Calculation of gammau & deltau
	double k9 = Xu > d ? h / (h - d) : std::min(h / (h - d), h / (d - Xu));
	double gammau = 2 * (exu + e2u) * k9 / tan(thetau * PI / 180.0);
	deltau = gammau * 2 * h;



	// YIELDING POINT

	// Calculation of Vy
	double fc1y = (0.2 + 0.3 * Xu / h) * 0.33 * sqrt(fc);
	double k10 = ((Ast * fyt) / (b * s) + fc1y) * b * dv;
	double k11 = (400 - thetau) / 360;
	double k12 = (45 - thetau) / (36 * vu);
	vy = (k10 * k12) > (k11 * k11 / 4) ? abs(k11 / (2 * k12)) : (k11 / 2 - sqrt(k11 * k11 / 4 - k10 * k12)) / k12;
	vy = std::min(vy,0.999 * vu);

	// Calculation of exy
	double exy_ten = (k1 * vy + k1p) / 750;
	double exy_comp = (k1 * vu + k1p) * Xu / (750 * (d - Xu));
	double exy = std::max(std::min(((exy_ten - exy_comp) / 2), exu), -1 * pRat * fc / (0.5 * Ec + Asl * Es / (b * d)));

	// Calculation of thetay
	double thetay = 45 - (45 - thetau) * (1.11 * vy / vu - 0.11);

	// Calculation of e2y
	double e2y = 0.9 * e2u * vy / vu;

	// Calculation of gammay & deltay
	double Xy = std::min(exy_comp / (exy_ten + exy_comp) * d, 0.99 * d);
	double k13 = std::min(h / (d - Xy), 1.2 * k9);
	double G = 0.4 * Ec;
	double gammay = std::max(std::min(2 * (exy + e2y) * k13 / tan(thetay * PI / 180), gammau), vy / (G * b * dv));
	gammay = std::min(gammay, gammau);
	deltay = gammay * 2 * h;


	// CRACKING POINT

	// Calculation of Vcr
	vcr = (0.33 * sqrt(fc) * b * h * h + pRat * fc * b * h * h) / (6 * alpha * dv);
	vcr = std::min(vcr, vy);

	// Calculation of gammacr & deltacr
	double gammacr = vcr / (G * b * dv);
	deltacr = gammacr * 2 * h;


	// Check deltay limits
	deltay = std::min(deltay, (vy-vcr) * (deltau-deltacr) / (vu - vcr) + deltacr);


	// FAILURE POINT

	// Calculation of Vf
	vf = 0;

	// Calculation of gammaf & deltaf
	double thetaf = pRat >= 0.2 ? 25 : 30;

	double gammaf = 0.04 * (1 + 1 / (tan(thetaf * PI / 180) * tan(thetaf * PI / 180))) *
		a  / ((1 / tan(thetaf * PI / 180) + pRat * tan(thetaf * PI / 180) * fc * b * h * s / (Ast * fyt * dv)) * 2 * h);
	gammaf = std::max(gammaf, gammau);
	deltaf = gammaf * 2 * h;


	// DEEP MEMBERS
	// Vud
	if (isDeep == 1)
	{
		double phiu = atan((h - Xu) / (a * 2));
		double cau = Xu - cover;
		double n = Es / Ec;
		double pz = Ast / (b * s);
		vu = vu * (1 + (cau * sin(phiu) * sin(phiu) * cos(phiu) * cos(phiu)) * (1 + n * pz / pow(sin(thetau * PI / 180), 4)) / (n * pz * dv / pow(tan(thetau * PI / 180), 2)));

		// Vyd
		double phiy = atan((h - Xy) / (a * 2));
		double cay = Xy - cover;
		vy = std::min(vy * (1 + (Ec * b * cay * sin(phiy) * sin(phiy) * cos(phiy) * cos(phiy)) * (gammay) / (vy)), 0.999 * vu);


		// deltad
		deltau = gammau * std::min(1.5 * a, std::min((a + d / 2), 2 * h));
		deltay = std::min(gammay * std::min(1.5 * a, std::min((a + d / 2), 2 * h)),deltau);
		deltacr = gammacr * std::min(1.5 * a, std::min((a + d / 2), 2 * h));
		deltaf = std::max(gammaf * std::min(1.5 * a, std::min((a + d / 2), 2 * h)),deltau);
	}

	vu = forcecf * vu;
	vy = forcecf * vy;
	vcr = forcecf * vcr;
	vf = forcecf * vf;
	deltau = dispcf * deltau;
	deltay = dispcf * deltay;
	deltacr = dispcf * deltacr;
	deltaf = dispcf * deltaf;

}


RCShearHinge::RCShearHinge() :UniaxialMaterial(0, MAT_TAG_RCShearHinge)
{
	tStrain = 0; tStress = 0; tTangent = 0;
	cStrain = 0; cStress = 0; cTangent = 0;
	vu = 0;
	deltau = 0;
	vy = 0;
	deltay = 0;
	vcr = 0;
	deltacr = 0;
	vf = 0;
	deltaf = 0;

}

RCShearHinge::~RCShearHinge()
{
	// does nothing
}

int
RCShearHinge::setTrialStrain(double s, double strainRate)
{
	//all variables to the last commit state
	tStrain = s;
	if (s <= deltacr && s >= -deltacr)
	{
		tTangent = vcr / deltacr;
		tStress = tTangent * s;
	}
	else if (s <= deltay && s >= -deltay)
	{
		tTangent = (vy - vcr) / (deltay - deltacr);
		if (s > 0)
		{
			tStress = vcr + tTangent * (s - deltacr);
		}
		else
		{
			tStress = -vcr - tTangent * (-1 * s - deltacr);
		}

	}
	else if (s <= deltau && s >= -deltau)
	{
		tTangent = (vu - vy) / (deltau - deltay);
		if (s > 0)
		{
			tStress = vy + tTangent * (s - deltay);
		}
		else
		{
			tStress = -vy - tTangent * (-1 * s - deltay);
		}
	}
	else if (s <= deltaf && s >= -deltaf) {
		tTangent = (vf - vu) / (deltaf - deltau);

		if (s > 0)
		{
			tStress = vu + tTangent * (s - deltau);
		}
		else
		{
			tStress = -vu - tTangent * (-1 * s - deltau);
		}
	}
	else // > deltaf
	{
		tStress = 0;
		tTangent = 1e-10;
	}

	return 0;
}

double RCShearHinge::getStress(void)
{

	return  (tStress);
}

double RCShearHinge::getTangent(void)
{
	return (tTangent);
}

double RCShearHinge::getInitialTangent(void)
{
	return (vcr / deltacr);
}


double RCShearHinge::getStrain(void)
{
	return tStrain;
}

double RCShearHinge::getStrainRate(void)
{
	return 0;
}

int RCShearHinge::commitState(void)
{
	cStrain = tStrain;
	cStress = tStress;
	cTangent = tTangent;
	return 0;
}

int RCShearHinge::revertToLastCommit(void)
{
	tStrain = cStrain;
	tStress = cStress;
	tTangent = cTangent;
	return 0;
}

int RCShearHinge::revertToStart(void)
{
	tStrain = tStress = 0;
	tTangent = getInitialTangent();
	return 0;
}

UniaxialMaterial*
RCShearHinge::getCopy(void)
{
	RCShearHinge* theCopy = new RCShearHinge();
	theCopy->cStrain = cStrain;
	theCopy->cStress = cStress;
	theCopy->cTangent = cTangent;
	theCopy->tStrain = tStrain;
	theCopy->tStress = tStress;
	theCopy->tTangent = tTangent;
	theCopy->vu = vu;
	theCopy->deltau = deltau;
	theCopy->vy = vy;
	theCopy->deltay = deltay;
	theCopy->vcr = vcr;
	theCopy->deltacr = deltacr;
	theCopy->vf = vf;
	theCopy->deltaf = deltaf;
	return theCopy;
}

int
RCShearHinge::sendSelf(int cTag, Channel& theChannel)
{
	int res = 0;
	static Vector data(15);		
	data(0) = this->getTag();
	data(1) = cStrain;
	data(2) = cStress;
	data(3) = cTangent;
	data(4) = tStrain;
	data(5) = tStress;
	data(6) = tTangent;
	data(7) = vu;
	data(8) = deltau;
	data(9) = vy;
	data(10) = deltay;
	data(11) = vcr;
	data(12) = deltacr;
	data(13) = vf;
	data(14) = deltaf;

	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "RCShearHinge::sendSelf() - failed to send data\n";

	return res;
}

int
RCShearHinge::recvSelf(int cTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int res = 0;
	static Vector data(15);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);

	if (res < 0) {
		opserr << "RCShearHinge::recvSelf() - failed to receive data\n";
		this->setTag(0);
	}
	else {
		this->setTag((int)data(0));
		cStrain = data(1);
		cStress = data(2);
		cTangent = data(3);
		tStrain = data(4);
		tStress = data(5);
		tTangent = data(6);
		vu = data(7);
		deltau = data(8);
		vy = data(9);
		deltay = data(10);
		vcr = data(11);
		deltacr = data(12);
		vf = data(13);
		deltaf = data(14);
	}

	return res;
}

void
RCShearHinge::Print(OPS_Stream& s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "RCShearHinge tag: " << this->getTag() << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"RCShearHinge\", ";
	}
}