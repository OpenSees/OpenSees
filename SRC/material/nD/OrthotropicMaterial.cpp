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
                                                                        
// $Revision: 1.4 $
// $Date: 2020-04-19 23:01:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/OrthotropicMaterial.h,v $

// Davide Raino, Massimo Petracca - ASDEA Software, Italy
//
// A Generic Orthotropic Material Wrapper that can convert any
// nonlinear isotropic material into an orthotropic one by means of tensor
// mapping
//

#include <OrthotropicMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>

void *OPS_OrthotropicMaterial(void)
{
	// check arguments
	int numArgs= OPS_GetNumRemainingInputArgs();
	if (numArgs < 17) {
		opserr << 
			"nDMaterial Orthotropic Error: Few arguments (< 17).\n"
			"nDMaterial Orthotropic $tag $theIsoMat $Ex $Ey $Ez $Gxy $Gyz $Gzx $vxy $vyz $vzx $Asigmaxx $Asigmayy $Asigmazz $Asigmaxyxy $Asigmayzyz $Asigmaxzxz.\n";
		return nullptr;
	}
	
	// get integer data
	int iData[2];
	int numData = 2;
	if (OPS_GetInt(&numData, iData) != 0)  {
		opserr << "nDMaterial Orthotropic Error: invalid nDMaterial tags.\n";
		return nullptr;
	}

	// get double data
	double dData[15];
	numData = 15;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "nDMaterial Orthotropic Error: invalid data for nDMaterial Orthotropic with tah " << iData[0] << ".\n";
		return nullptr;
	}

	// get the isotropic material to map
	NDMaterial *theIsoMaterial = OPS_getNDMaterial(iData[1]);
	if (theIsoMaterial == 0) {
		opserr << "WARNING: nDMaterial does not exist.\n";
		opserr << "nDMaterial: " << iData[1] << "\n";
		opserr << "nDMaterial Orthotropic: " << iData[0] << "\n";
		return nullptr;
	}

	// create the orthotropic wrapper
	NDMaterial* theOrthotropicMaterial = new OrthotropicMaterial(
		iData[0], 
		*theIsoMaterial,
		dData[0], dData[1], dData[2], dData[3],
		dData[4], dData[5], dData[6], dData[7],
		dData[8], dData[9], dData[10], dData[11],
		dData[12], dData[13], dData[14]);
	if (theOrthotropicMaterial == 0) {
		opserr << "nDMaterial Orthotropic Error: failed to allocate a new material.\n";
		return nullptr;
	}

	// done
	return theOrthotropicMaterial;
}

OrthotropicMaterial::OrthotropicMaterial(
	int tag, 
	NDMaterial &theIsoMat,
	double Ex, double Ey, double Ez, double Gxy, double Gyz, double Gzx,
	double vxy, double vyz, double vzx,
	double Asigmaxx, double Asigmayy, double Asigmazz, double Asigmaxyxy, double Asigmayzyz, double Asigmaxzxz)
	: NDMaterial(tag, ND_TAG_OrthotropicMaterial)
{
	// copy the isotropic material
	theIsotropicMaterial = theIsoMat.getCopy("ThreeDimensional");
	if (theIsotropicMaterial == 0) {
		opserr << "nDMaterial Orthotropic Error: failed to get a (3D) copy of the isotropic material\n";
		exit(-1);
	}

	// compute the initial orthotropic constitutive tensor
	static Matrix C0(6, 6);
	C0.Zero();
	double vyx = vxy * Ey / Ex;
	double vzy = vyz * Ez / Ey;
	double vxz = vzx * Ex / Ez;
	double d = (1.0 - vxy * vyx - vyz * vzy - vzx * vxz - 2.0*vxy*vyz*vzx) / (Ex*Ey*Ez);
	C0(0, 0) = (1.0 - vyz * vzy) / (Ey*Ez*d);
	C0(1, 1) = (1.0 - vzx * vxz) / (Ez*Ex*d);
	C0(2, 2) = (1.0 - vxy * vyx) / (Ex*Ey*d);
	C0(1, 0) = (vxy + vxz * vzy) / (Ez*Ex*d);
	C0(0, 1) = C0(1, 0);
	C0(2, 0) = (vxz + vxy * vyz) / (Ex*Ey*d);
	C0(0, 2) = C0(2, 0);
	C0(2, 1) = (vyz + vxz * vyx) / (Ex*Ey*d);
	C0(1, 2) = C0(2, 1);
	C0(3, 3) = Gxy;
	C0(4, 4) = Gyz;
	C0(5, 5) = Gzx;

	// compute the Asigma and its inverse
	if (Asigmaxx <= 0 || Asigmayy <= 0 || Asigmazz <= 0 || Asigmaxyxy <= 0 || Asigmayzyz <= 0 || Asigmaxzxz <= 0) {
		opserr << "nDMaterial Orthotropic Error: Asigma11, Asigma22, Asigma33, Asigma12, Asigma23, Asigma13 must be greater than 0.\n";
		exit(-1);
	}
	static Matrix Asigma(6, 6);
	Asigma.Zero();
	Asigma(0, 0) = Asigmaxx;
	Asigma(1, 1) = Asigmayy;
	Asigma(2, 2) = Asigmazz;
	Asigma(3, 3) = Asigmaxyxy;
	Asigma(4, 4) = Asigmayzyz;
	Asigma(5, 5) = Asigmaxzxz;
	for (int i = 0; i < 6; ++i)
		Asigma_inv(i) = 1.0 / Asigma(i, i);

	// coompute the initial isotropic constitutive tensor and its inverse
	static Matrix C0iso(6, 6);
	static Matrix C0iso_inv(6, 6);
	C0iso = theIsotropicMaterial->getInitialTangent();
	int res = C0iso.Invert(C0iso_inv);
	if (res < 0) {
		opserr << "nDMaterial Orthotropic Error: the isotropic material gave a singular initial tangent.\n";
		exit(-1);
	}

	// compute the strain tensor map inv(C0_iso) * Asigma * C0_ortho
	static Matrix Asigma_C0(6, 6);
	Asigma_C0.addMatrixProduct(0.0, Asigma, C0, 1.0);
	Aepsilon.addMatrixProduct(0.0, C0iso_inv, Asigma_C0, 1.0);
}

OrthotropicMaterial::OrthotropicMaterial()
	: NDMaterial(0, ND_TAG_OrthotropicMaterial)
{
}

OrthotropicMaterial::~OrthotropicMaterial()
{ 
	if (theIsotropicMaterial)
		delete theIsotropicMaterial;
}

double OrthotropicMaterial::getRho(void)
{
	return theIsotropicMaterial->getRho();
}

int OrthotropicMaterial::setTrialStrain(const Vector & strain)
{
	// strain in orthotropic space
	epsilon = strain;

	// move to isotropic space
	static Vector eps_iso(6);
	eps_iso.addMatrixVector(0.0, Aepsilon, epsilon, 1.0);

	// call isotropic material
	int res = theIsotropicMaterial->setTrialStrain(eps_iso);
	if (res != 0) {
		opserr << "nDMaterial Orthotropic Error: the isotropic material failed in setTrialStrain.\n";
		return res;
	}
	return 0;
}

const Vector &OrthotropicMaterial::getStrain(void)
{
	return epsilon;
}

const Vector &OrthotropicMaterial::getStress(void)
{
	// stress in isotropic space
	const Vector& sigma_iso = theIsotropicMaterial->getStress();

	// move to orthotropic space
	static Vector sigma(6);
	for (int i = 0; i < 6; ++i)
		sigma(i) = Asigma_inv(i) * sigma_iso(i);
	return sigma;
}

const Matrix &OrthotropicMaterial::getTangent(void)
{
	// tensor in isotropic space
	const Matrix &C_iso = theIsotropicMaterial->getTangent();

	// compute orthotripic tangent
	static Matrix C(6, 6);
	static Matrix temp(6, 6);
	static Matrix invAsigma(6, 6);
	invAsigma.Zero();
	for (int i = 0; i < 6; ++i)
		invAsigma(i, i) = Asigma_inv(i);
	temp.addMatrixProduct(0.0, C_iso, Aepsilon, 1.0);
	C.addMatrixProduct(0.0, invAsigma, temp, 1.0);
	return C;
}

const Matrix &OrthotropicMaterial::getInitialTangent(void)
{
	// tensor in isotropic space
	const Matrix& C_iso = theIsotropicMaterial->getInitialTangent();

	// compute orthotripic tangent
	static Matrix C(6, 6);
	static Matrix temp(6, 6);
	static Matrix invAsigma(6, 6);
	invAsigma.Zero();
	for (int i = 0; i < 6; ++i)
		invAsigma(i, i) = Asigma_inv(i);
	temp.addMatrixProduct(0.0, C_iso, Aepsilon, 1.0);
	C.addMatrixProduct(0.0, invAsigma, temp, 1.0);
	return C;;
}

int OrthotropicMaterial::commitState(void)
{
	return theIsotropicMaterial->commitState();
}

int OrthotropicMaterial::revertToLastCommit(void)
{
	return theIsotropicMaterial->revertToLastCommit();
}

int OrthotropicMaterial::revertToStart(void)
{
	return theIsotropicMaterial->revertToStart();
}

NDMaterial * OrthotropicMaterial::getCopy(void)
{
	OrthotropicMaterial *theCopy = new OrthotropicMaterial();
	theCopy->setTag(getTag());
	theCopy->theIsotropicMaterial = theIsotropicMaterial->getCopy("ThreeDimensional");
	theCopy->epsilon = epsilon;
	theCopy->Aepsilon = Aepsilon;
	theCopy->Asigma_inv = Asigma_inv;
	return theCopy;
}

NDMaterial* OrthotropicMaterial::getCopy(const char* code)
{
	if (strcmp(code, "ThreeDimensional") == 0)
		return getCopy();
	return NDMaterial::getCopy(code);
}

const char* OrthotropicMaterial::getType(void) const
{
	return "ThreeDimensional";
}

int OrthotropicMaterial::getOrder(void) const
{
	return 6;
}

void OrthotropicMaterial::Print(OPS_Stream &s, int flag)
{
	s << "Orthotropic Material, tag: " << this->getTag() << "\n";
}

int OrthotropicMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	// result
	int res = 0;

	// data
	static Vector data(48);
	int counter = 0;
	// store int values
	data(counter++) = static_cast<double>(getTag());
	data(counter++) = static_cast<double>(theIsotropicMaterial->getClassTag());
	int matDbTag = theIsotropicMaterial->getDbTag();
	if (matDbTag == 0) {
		matDbTag = theChannel.getDbTag();
		theIsotropicMaterial->setDbTag(matDbTag);
	}
	data(counter++) = static_cast<double>(matDbTag);
	// store internal variables
	for (int i = 0; i < 6; ++i)
		data(counter++) = epsilon(i);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			data(counter++) = Aepsilon(i, j);
	for (int i = 0; i < 6; ++i)
		data(counter++) = Asigma_inv(i);
	// send data
	res = theChannel.sendVector(getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "nDMaterial Orthotropic Error: failed to send vector data\n";
		return res;
	}

	// now send the materials data
	res = theIsotropicMaterial->sendSelf(commitTag, theChannel);
	if (res < 0) {
		opserr << "nDMaterial Orthotropic Error: failed to send the isotropic material\n";
		return res;
	}

	// done
	return res;
}

int OrthotropicMaterial::recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
{
	// result
	int res = 0;

	// data
	static Vector data(48);
	int counter = 0;

	// receive data
	res = theChannel.recvVector(getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "nDMaterial Orthotropic Error: failed to send vector data\n";
		return res;
	}
	// get int values
	setTag(static_cast<int>(data(counter++)));
	int matClassTag = static_cast<int>(data(counter++));
	// if the associated material has not yet been created or is of the wrong type
	// create a new material for recvSelf later
	if ((theIsotropicMaterial == nullptr) || (theIsotropicMaterial->getClassTag() != matClassTag)) {
		if (theIsotropicMaterial)
			delete theIsotropicMaterial;
		theIsotropicMaterial = theBroker.getNewNDMaterial(matClassTag);
		if (theIsotropicMaterial == nullptr) {
			opserr << "nDMaterial Orthotropic Error: failed to get a material of type: " << matClassTag << endln;
			return -1;
		}
	}
	theIsotropicMaterial->setDbTag(static_cast<int>(data(counter++)));
	// store internal variables
	for (int i = 0; i < 6; ++i)
		epsilon(i) = data(counter++);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			Aepsilon(i, j) = data(counter++);
	for (int i = 0; i < 6; ++i)
		Asigma_inv(i) = data(counter++);

	// now receive the associated materials data
	res = theIsotropicMaterial->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		opserr << "nDMaterial Orthotropic Error: failed to receive the isotropic material\n";
		return res;
	}

	// done
	return res;
}

int OrthotropicMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	// forward to the adapted (isotropic) material
	return theIsotropicMaterial->setParameter(argv, argc, param);
}

Response* OrthotropicMaterial::setResponse(const char** argv, int argc, OPS_Stream& s)
{
	if (argc > 0) {
		if (strcmp(argv[0], "stress") == 0 || 
			strcmp(argv[0], "stresses") == 0 || 
			strcmp(argv[0], "strain") == 0 || 
			strcmp(argv[0], "strains") == 0 ||
			strcmp(argv[0], "Tangent") == 0 || 
			strcmp(argv[0], "tangent") == 0) {
			// stresses, strain and tangent should be those of this adapter (orthotropic)
			return NDMaterial::setResponse(argv, argc, s);
		}
		else {
			// any other response should be obtained from the adapted (isotropic) material
			return theIsotropicMaterial->setResponse(argv, argc, s);
		}
	}
	return NDMaterial::setResponse(argv, argc, s);
}
