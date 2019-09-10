#include "stdafx.h"
#include "ZeroLengthElementWrapper.h"

using namespace System::Runtime::InteropServices;
using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;



//CoupledZeroLengthWrapper::CoupledZeroLengthWrapper(int tag,
//	int Nd1, int Nd2,
//	UniaxialMaterialWrapper^ theMaterial,
//	int direction1, int direction2,
//	int doRayleighDamping) {
//
//	_Element = new CoupledZeroLength( tag,
//		 Nd1,  Nd2,
//		*theMaterial->_UniaxialMaterial,
//		 direction1,  direction2,
//		 doRayleighDamping);
//}

ZeroLengthWrapper::ZeroLengthWrapper(int tag,
	int dimension,
	int Nd1, int Nd2,
	VectorWrapper^ x,
	VectorWrapper^ yprime,
	UniaxialMaterialWrapper^ theMaterial,
	int direction,
	int doRayleighDamping) {

	_Element = new ZeroLength(tag,
		dimension,
		Nd1, Nd2,
		*x->_Vector,
		*yprime->_Vector,
		*theMaterial->_UniaxialMaterial,
		direction,
		doRayleighDamping);
}

ZeroLengthWrapper::ZeroLengthWrapper(int tag,
	int dimension,
	int Nd1, int Nd2,
	VectorWrapper^ x,
	VectorWrapper^ yprime,
	UniaxialMaterialWrapper^ theMaterial,
	UniaxialMaterialWrapper^ theDampingMaterial,
	int direction) {

	_Element = new ZeroLength(tag,
		dimension,
		Nd1, Nd2,
		*x->_Vector,
		*yprime->_Vector,
		*theMaterial->_UniaxialMaterial,
		*theDampingMaterial->_UniaxialMaterial,
		direction);
}

UniaxialMaterial** array2pointer3(array<UniaxialMaterialWrapper^>^ theMaterials) {
	UniaxialMaterial** mats = new UniaxialMaterial*[theMaterials->Length];
	for (int i = 0; i < theMaterials->Length; i++)
	{
		mats[i] = theMaterials[i]->_UniaxialMaterial;
	}

	return mats;
}

ZeroLengthWrapper::ZeroLengthWrapper(int tag,
	int dimension,
	int Nd1, int Nd2,
	VectorWrapper^ x,
	VectorWrapper^ yprime,
	int n1dMat,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	IDWrapper^ direction,
	int doRaylieghDamping) {

	_Element = new ZeroLength(tag,
		dimension,
		Nd1, Nd2,
		*x->_Vector,
		*yprime->_Vector,
		n1dMat,
		array2pointer3(theMaterials),
		*direction->_ID,
		doRaylieghDamping);
}

ZeroLengthWrapper::ZeroLengthWrapper(int tag,
	int dimension,
	int Nd1, int Nd2,
	VectorWrapper^ x,
	VectorWrapper^ yprime,
	int n1dMat,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	array<UniaxialMaterialWrapper^>^ theDampMaterial,
	IDWrapper^ direction,
	int doRaylieghDamping) {

	_Element = new ZeroLength(tag,
		dimension,
		Nd1, Nd2,
		*x->_Vector,
		*yprime->_Vector,
		n1dMat,
		array2pointer3(theMaterials),
		array2pointer3(theDampMaterial),
		*direction->_ID,
		doRaylieghDamping);
}

ZeroLengthContact2DWrapper::ZeroLengthContact2DWrapper(int tag, int Nd1, int Nd2,
	double Kn, double Kt, double fRatio,
	VectorWrapper^ normal) {

	_Element = new ZeroLengthContact2D(tag, Nd1, Nd2,
		Kn, Kt, fRatio,
		*normal->_Vector);
}

ZeroLengthContact3DWrapper::ZeroLengthContact3DWrapper(int tag,
	int Nd1, int Nd2,
	int direction, double Kn, double Kt, double fRatio, double c,
	double originX, double originY) {

	_Element = new ZeroLengthContact3D(tag,
		Nd1, Nd2,
		direction, Kn, Kt, fRatio, c,
		originX, originY);
}

ZeroLengthContactNTS2DWrapper::ZeroLengthContactNTS2DWrapper(int tag, int sNdNum, int mNdNum, IDWrapper^ Nodes,
	double Kn, double Kt, double fRatio) {

	_Element = new ZeroLengthContactNTS2D(tag, sNdNum, mNdNum, *Nodes->_ID,
		Kn, Kt, fRatio);
}

ZeroLengthImpact3DWrapper::ZeroLengthImpact3DWrapper(int tag,
	int Nd1, int Nd2,
	int direction,
	double initGapInput, double fRatio, double Kt,
	double Kn, double Kn2Input, double Delta_yInput,
	double c) {

	_Element = new ZeroLengthImpact3D(tag,
		Nd1, Nd2,
		direction,
		initGapInput, fRatio, Kt,
		Kn, Kn2Input, Delta_yInput,
		c);
}

ZeroLengthInterface2DWrapper::ZeroLengthInterface2DWrapper(int tag, int sNdNum, int mNdNum, int sDof, int mDof, IDWrapper^ Nodes,
	double Kn, double Kt, double fRatio) {

	_Element = new ZeroLengthInterface2D(tag, sNdNum, mNdNum, sDof, mDof, *Nodes->_ID,
		Kn, Kt, fRatio);
}

ZeroLengthNDWrapper::ZeroLengthNDWrapper(int tag,
int dimension,
int Nd1, int Nd2,
VectorWrapper^ x,
VectorWrapper^ yprime,
NDMaterialWrapper^ theNDMaterial) {

	_Element = new ZeroLengthND( tag,
		 dimension,
		 Nd1,  Nd2,
		*x->_Vector,
		*yprime->_Vector,
		*theNDMaterial->_NDMaterial);
}

ZeroLengthNDWrapper::ZeroLengthNDWrapper(int tag,
	int dimension,
	int Nd1, int Nd2,
	VectorWrapper^ x,
	VectorWrapper^ yprime,
	NDMaterialWrapper^ theNDMaterial, UniaxialMaterialWrapper^ the1DMaterial) {

	_Element = new ZeroLengthND( tag,
		 dimension,
		 Nd1,  Nd2,
		*x->_Vector,
		*yprime->_Vector,
		*theNDMaterial->_NDMaterial, *the1DMaterial->_UniaxialMaterial);
}

ZeroLengthRockingWrapper::ZeroLengthRockingWrapper(int tag,
	int dimension,
	int Nd1, int Nd2,
	VectorWrapper^ x,
	VectorWrapper^ yprime,
	double kr, double R,
	double theta, double kappa,
	double xi, double dispTol, double velTol) {

	_Element = new ZeroLengthRocking(tag,
		dimension,
		Nd1, Nd2,
		*x->_Vector,
		*yprime->_Vector,
		kr, R,
		theta, kappa,
		xi, dispTol, velTol);
}

ZeroLengthSectionWrapper::ZeroLengthSectionWrapper(int tag,
	int dimension,
	int Nd1, int Nd2,
	VectorWrapper^ x,
	VectorWrapper^ yprime,
	SectionForceDeformationWrapper^ theSection,
	int doRayleighDamping) {

	_Element = new ZeroLengthSection(tag,
		dimension,
		Nd1, Nd2,
		*x->_Vector,
		*yprime->_Vector,
		*theSection->_SectionForceDeformation,
		doRayleighDamping);
}





