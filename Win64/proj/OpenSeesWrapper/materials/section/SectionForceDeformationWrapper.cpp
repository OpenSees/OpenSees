#include "stdafx.h"
#include "SectionForceDeformationWrapper.h"

using namespace System;
using namespace System::Collections::Generic;
using namespace OpenSees;
using namespace OpenSees::Materials;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::NDMaterials;

double* array2pointer6(array<double>^ arr) {
	double* rets = new double[arr->Length];
	for (int i = 0; i < arr->Length; i++)
	{
		rets[i] = arr[i];
	}

	return rets;
}

UniaxialMaterial** array2pointer6(array<UniaxialMaterialWrapper^>^ theMaterials) {
	UniaxialMaterial** mats = new UniaxialMaterial*[theMaterials->Length];
	for (int i = 0; i < theMaterials->Length; i++)
	{
		mats[i] = theMaterials[i]->_UniaxialMaterial;
	}

	return mats;
}

NDMaterial** array2pointer6(array<NDMaterialWrapper^>^ theMaterials) {
	NDMaterial** mats = new NDMaterial*[theMaterials->Length];
	for (int i = 0; i < theMaterials->Length; i++)
	{
		mats[i] = theMaterials[i]->_NDMaterial;
	}

	return mats;
}

SectionForceDeformation** array2pointer6(array<SectionForceDeformationWrapper^>^ theMaterials) {
	SectionForceDeformation** mats = new SectionForceDeformation*[theMaterials->Length];
	for (int i = 0; i < theMaterials->Length; i++)
	{
		mats[i] = theMaterials[i]->_SectionForceDeformation;
	}

	return mats;
}

Fiber** array2pointer6(array<FiberWrapper^>^ fibers) {
	Fiber** mats = new Fiber*[fibers->Length];
	for (int i = 0; i < fibers->Length; i++)
	{
		mats[i] = fibers[i]->_Fiber;
	}
	return mats;
}

SectionForceDeformationWrapper::SectionForceDeformationWrapper()
{

}

FiberSection3dWrapper::FiberSection3dWrapper(int tag, array<FiberWrapper^>^ fibers, UniaxialMaterialWrapper^ GJ)
{
	//QuadPatchFiber3d **_QuadPatchFiber3ds = 0;
	if (GJ != nullptr)
		_GJ = GJ->_UniaxialMaterial;

	_SectionForceDeformation = new FiberSection3d(tag, fibers->Length, array2pointer6(fibers), _GJ);
}

FiberSection3dWrapper::FiberSection3dWrapper(int tag, FiberSectionReprWrapper ^ fiberSectionRepr,
	Dictionary<int, UniaxialMaterialWrapper^>^ materials, array<UniaxialFiber3dWrapper^>^ fibers, UniaxialMaterialWrapper^ GJ)
{
	int numPatches = fiberSectionRepr->_FiberSectionRepr->getNumPatches();
	Patch** _patches = fiberSectionRepr->_FiberSectionRepr->getPatches();
	int numReinfLayers = fiberSectionRepr->_FiberSectionRepr->getNumReinfLayers();
	ReinfLayer** _reinfLayers = fiberSectionRepr->_FiberSectionRepr->getReinfLayers();
	int numFibers = 0;
	for (int i = 0; i < numPatches; i++)
		numFibers += _patches[i]->getNumCells();
	for (int i = 0; i < numReinfLayers; i++)
		numFibers += _reinfLayers[i]->getNumReinfBars();

	int numAddFibers = 0;
	if (fibers != nullptr)
		numAddFibers = fibers->Length;
	numFibers += numAddFibers;

	Fiber** _fibers = new Fiber*[numFibers];
	int currentFiber = 0;

	for (int i = 0; i < numPatches; i++)
	{
		int numCells = _patches[i]->getNumCells();
		Cell** _cells = _patches[i]->getCells();
		int theMatTag = _patches[i]->getMaterialID();

		if (!materials->ContainsKey(theMatTag))
			throw gcnew System::Exception("UniaxialMaterial with tag" + theMatTag.ToString() + " not found ?!");

		UniaxialMaterial* theMat = materials[theMatTag]->_UniaxialMaterial;

		for (int j = 0; j < numCells; j++)
		{
			_fibers[currentFiber] = new UniaxialFiber3d(currentFiber, *theMat, _cells[j]->getArea(),
				_cells[j]->getCentroidPosition());
			currentFiber++;
		}
	}
	for (int i = 0; i < numReinfLayers; i++)
	{
		int numReinfBars = _reinfLayers[i]->getNumReinfBars();
		ReinfBar* _reinfBars = _reinfLayers[i]->getReinfBars();
		int theMatTag = _reinfLayers[i]->getMaterialID();

		if (!materials->ContainsKey(theMatTag))
			throw gcnew System::Exception("UniaxialMaterial with tag" + theMatTag.ToString() + " not found ?!");

		UniaxialMaterial* theMat = materials[theMatTag]->_UniaxialMaterial;

		for (int j = 0; j < numReinfBars; j++)
		{
			_fibers[currentFiber] = new UniaxialFiber3d(currentFiber, *theMat, _reinfBars[i].getArea(),
				_reinfBars[i].getPosition());
			currentFiber++;
		}
	}
	for (int i = 0; i < numAddFibers; i++)
	{
		_fibers[currentFiber] = fibers[i]->_Fiber;
		currentFiber++;
	}


	_GJ = GJ->_UniaxialMaterial;
	_SectionForceDeformation = new FiberSection3d(tag, numFibers, _fibers, _GJ);
}


FiberSection2dWrapper::FiberSection2dWrapper(int tag, array<FiberWrapper^>^ fibers)
{
	_SectionForceDeformation = new FiberSection2d(tag, fibers->Length, array2pointer6(fibers));

}

FiberSection2dWrapper::FiberSection2dWrapper(int tag, FiberSectionReprWrapper ^ fiberSectionRepr,
	Dictionary<int, UniaxialMaterialWrapper^>^ materials, array<UniaxialFiber2dWrapper^>^ fibers)
{
	int numPatches = fiberSectionRepr->_FiberSectionRepr->getNumPatches();
	Patch** _patches = fiberSectionRepr->_FiberSectionRepr->getPatches();
	int numReinfLayers = fiberSectionRepr->_FiberSectionRepr->getNumReinfLayers();
	ReinfLayer** _reinfLayers = fiberSectionRepr->_FiberSectionRepr->getReinfLayers();
	int numFibers = 0;
	for (int i = 0; i < numPatches; i++)
		numFibers += _patches[i]->getNumCells();
	for (int i = 0; i < numReinfLayers; i++)
		numFibers += _reinfLayers[i]->getNumReinfBars();
	int numAddFibers = 0;
	if (fibers != nullptr)
		numAddFibers = fibers->Length;
	numFibers += numAddFibers;

	Fiber** _fibers = new Fiber*[numFibers];
	int currentFiber = 0;

	for (int i = 0; i < numPatches; i++)
	{
		int numCells = _patches[i]->getNumCells();
		Cell** _cells = _patches[i]->getCells();
		int theMatTag = _patches[i]->getMaterialID();

		if (!materials->ContainsKey(theMatTag))
			throw gcnew System::Exception("UniaxialMaterial with tag" + theMatTag.ToString() + " not found ?!");

		UniaxialMaterial* theMat = materials[theMatTag]->_UniaxialMaterial;

		for (int j = 0; j < numCells; j++)
		{
			_fibers[currentFiber] = new UniaxialFiber2d(currentFiber, *theMat, _cells[j]->getArea(),
				_cells[j]->getCentroidPosition()[0]);
			currentFiber++;
		}
	}

	for (int i = 0; i < numReinfLayers; i++)
	{
		int numReinfBars = _reinfLayers[i]->getNumReinfBars();
		ReinfBar* _reinfBars = _reinfLayers[i]->getReinfBars();
		int theMatTag = _reinfLayers[i]->getMaterialID();

		if (!materials->ContainsKey(theMatTag))
			throw gcnew System::Exception("UniaxialMaterial with tag" + theMatTag.ToString() + " not found ?!");

		UniaxialMaterial* theMat = materials[theMatTag]->_UniaxialMaterial;

		for (int j = 0; j < numReinfBars; j++)
		{
			_fibers[currentFiber] = new UniaxialFiber2d(currentFiber, *theMat, _reinfBars[i].getArea(),
				_reinfBars[i].getPosition()[0]);
			currentFiber++;
		}
	}

	for (int i = 0; i < numAddFibers; i++)
	{
		_fibers[currentFiber] = fibers[i]->_Fiber;
		currentFiber++;
	}

	_SectionForceDeformation = new FiberSection2d(tag, numFibers, _fibers);
}

ElasticSection3dWrapper::ElasticSection3dWrapper(int tag, double E, double A, double Iz,
	double Iy, double G, double J)
{
	_SectionForceDeformation = new ElasticSection3d(tag, E, A, Iz, Iy, G, J);
}

FiberSection2dIntWrapper::FiberSection2dIntWrapper(int tag,
	array<FiberWrapper^>^ fibers,
	array<FiberWrapper^>^ Hfibers,
	int NStrip1,
	double tavg1,
	int NStrip2,
	double tavg2,
	int NStrip3,
	double tavg3)
{
	int numFibers = fibers->Length;
	int numHFibers = Hfibers->Length;

	Fiber** _fibers = new Fiber*[numFibers];
	Fiber** _Hfibers = new Fiber*[numHFibers];

	for (int i = 0; i < numFibers; i++)
		_fibers[i] = fibers[i]->_Fiber;
	for (int i = 0; i < numHFibers; i++)
		_Hfibers[i] = Hfibers[i]->_Fiber;

	_SectionForceDeformation = new FiberSection2dInt(tag,
		numFibers,
		_fibers,
		numHFibers,
		_Hfibers,
		NStrip1,
		tavg1,
		NStrip2,
		tavg2,
		NStrip3,
		tavg3);
}

YS_Section2D01Wrapper::YS_Section2D01Wrapper(int tag, double E, double A, double I,
	YieldSurface_BCWrapper^ ptrys,
	bool use_kr) {
	_SectionForceDeformation = new YS_Section2D01(tag, E, A, I, ptrys->_YieldSurface_BC, use_kr);
}

YS_Section2D02Wrapper::YS_Section2D02Wrapper(int tag,
	double E, double A, double I,
	double theta_p_max,
	YieldSurface_BCWrapper^ ptrys,
	bool use_kr) {
	_SectionForceDeformation = new YS_Section2D02(tag, E, A, I, theta_p_max, ptrys->_YieldSurface_BC, use_kr);
}

BidirectionalWrapper::BidirectionalWrapper(int tag, double E, double sigY, double Hiso, double Hkin,
	int code1,
	int code2) {
	_SectionForceDeformation = new Bidirectional(tag, E, sigY, Hiso, Hkin, code1, code2);
}

ElasticMembranePlateSectionWrapper::ElasticMembranePlateSectionWrapper(int    tag,
	double E,
	double nu,
	double h,
	double rho) {
	_SectionForceDeformation = new ElasticMembranePlateSection(tag, E, nu, h, rho);
}

ElasticPlateSectionWrapper::ElasticPlateSectionWrapper(int    tag,
	double E,
	double nu,
	double h) {
	_SectionForceDeformation = new ElasticPlateSection(tag, E, nu, h);
}

ElasticSection2dWrapper::ElasticSection2dWrapper(int tag, double E, double A, double I) {
	_SectionForceDeformation = new ElasticSection2d(tag, E, A, I);
}



ElasticShearSection2dWrapper::ElasticShearSection2dWrapper(int tag, double E, double A, double I,
	double G, double alpha) {
	_SectionForceDeformation = new ElasticShearSection2d(tag, E, A, I, G, alpha);
}

ElasticShearSection3dWrapper::ElasticShearSection3dWrapper(int tag, double E, double A, double Iz,
	double Iy, double G, double J, double alphaY, double alphaZ) {
	_SectionForceDeformation = new ElasticShearSection3d(tag, E, A, Iz, Iy, G, J, alphaY, alphaZ);
}

ElasticTubeSection3dWrapper::ElasticTubeSection3dWrapper(int tag, double E, double d, double tw, double G) {
	_SectionForceDeformation = new ElasticTubeSection3d(tag, E, d, tw, G);
}

ElasticWarpingShearSection2dWrapper::ElasticWarpingShearSection2dWrapper(int tag, double E, double A, double I,
	double G, double alpha, double J, double B, double C) {
	_SectionForceDeformation = new ElasticWarpingShearSection2d(tag, E, A, I, G, alpha, J, B, C);
}

Elliptical2Wrapper::Elliptical2Wrapper(int tag, double E1, double E2, double sigY1, double sigY2,
	double Hiso, double Hkin1, double Hkin2,
	int c1, int c2) {
	_SectionForceDeformation = new Elliptical2(tag, E1, E2, sigY1, sigY2, Hiso, Hkin1, Hkin2, c1, c2);
}

FiberSection2dThermalWrapper::FiberSection2dThermalWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers) {
	_SectionForceDeformation = new FiberSection2dThermal(tag, fibers->Length, array2pointer6(fibers));
}

FiberSection3dThermalWrapper::FiberSection3dThermalWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers) {
	_SectionForceDeformation = new FiberSection3dThermal(tag, fibers->Length, array2pointer6(fibers));
}

FiberSectionGJWrapper::FiberSectionGJWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers, double GJ) {
	_SectionForceDeformation = new FiberSectionGJ(tag, fibers->Length, array2pointer6(fibers), GJ);
}

FiberSectionGJThermalWrapper::FiberSectionGJThermalWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers, double GJ) {
	_SectionForceDeformation = new FiberSectionGJThermal(tag, fibers->Length, array2pointer6(fibers), GJ);
}

GenericSection1dWrapper::GenericSection1dWrapper(int tag, UniaxialMaterialWrapper^ m, int code) {
	_SectionForceDeformation = new GenericSection1d(tag, *m->_UniaxialMaterial, code);
}

Isolator2springWrapper::Isolator2springWrapper(int tag, double tol_in, double k1_in, double Fy_in, double kb_in, double kvo_in,
	double hb_in, double Pe_in, double po_in) {
	_SectionForceDeformation = new Isolator2spring(tag, tol_in, k1_in, Fy_in, kb_in, kvo_in, hb_in, Pe_in, po_in);
}

LayeredShellFiberSectionWrapper::LayeredShellFiberSectionWrapper(int tag,
	int iLayers,
	array<double>^ thickness,
	array<NDMaterialWrapper^>^ fibers) {
	_SectionForceDeformation = new LayeredShellFiberSection(tag, iLayers, array2pointer6(thickness), array2pointer6(fibers));
}

LayeredShellFiberSectionThermalWrapper::LayeredShellFiberSectionThermalWrapper(int tag,
	int iLayers,
	array<double>^ thickness,
	array<NDMaterialWrapper^>^ fibers) {
	_SectionForceDeformation = new LayeredShellFiberSectionThermal(tag, iLayers, array2pointer6(thickness), array2pointer6(fibers));
}

MembranePlateFiberSectionWrapper::MembranePlateFiberSectionWrapper(int tag,
	double thickness,
	NDMaterialWrapper^ fibers) {
	_SectionForceDeformation = new MembranePlateFiberSection(tag, thickness, *fibers->_NDMaterial);
}

MembranePlateFiberSectionThermalWrapper::MembranePlateFiberSectionThermalWrapper(int tag,
	double thickness,
	NDMaterialWrapper^ fibers) {
	_SectionForceDeformation = new MembranePlateFiberSectionThermal(tag, thickness, *fibers->_NDMaterial);
}

NDFiberSection2dWrapper::NDFiberSection2dWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers, double a) {
	_SectionForceDeformation = new NDFiberSection2d(tag, fibers->Length, array2pointer6(fibers), a);
}

NDFiberSection3dWrapper::NDFiberSection3dWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers, double a) {
	_SectionForceDeformation = new NDFiberSection3d(tag, numFibers, array2pointer6(fibers), a);
}

NDFiberSectionWarping2dWrapper::NDFiberSectionWarping2dWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers, double a) {
	_SectionForceDeformation = new NDFiberSectionWarping2d(tag, numFibers, array2pointer6(fibers), a);
}

ParallelSectionWrapper::ParallelSectionWrapper(int tag, int numSections, array<SectionForceDeformationWrapper^>^ theSections) {
	_SectionForceDeformation = new ParallelSection(tag, numSections, array2pointer6(theSections));
}

SectionAggregatorWrapper::SectionAggregatorWrapper(int tag, SectionForceDeformationWrapper^ theSection,
	int numAdditions, array<UniaxialMaterialWrapper^>^ theAdditions,
	IDWrapper^ code) {
	_SectionForceDeformation = new SectionAggregator(tag, *theSection->_SectionForceDeformation, theAdditions->Length, array2pointer6(theAdditions), *code->_ID);
}

SectionAggregatorWrapper::SectionAggregatorWrapper(int tag, int numAdditions, array<UniaxialMaterialWrapper^>^ theAdditions,
	IDWrapper^ code) {
	_SectionForceDeformation = new SectionAggregator(tag, theAdditions->Length, array2pointer6(theAdditions), *code->_ID);
}


SectionAggregatorWrapper::SectionAggregatorWrapper(int tag, SectionForceDeformationWrapper^ theSection, UniaxialMaterialWrapper^ theAdditions,
	int code) {
	_SectionForceDeformation = new

		SectionAggregator(tag, *theSection->_SectionForceDeformation, *theAdditions->_UniaxialMaterial, code);
}

//FiberSection3dWrapper::FiberSection3dWrapper(int tag, array<FiberWrapper^>^ fibers, array<QuadPatchFiber3dWrapper^>^ quadPatchFiber3ds, UniaxialMaterialWrapper^ GJ)
//{
//	int numFibers = fibers->Length;
//	Fiber** _fibers = new Fiber*[numFibers];
//	for (int i = 0; i < numFibers; i++)
//		_fibers[i] = fibers[i]->_Fiber;
//
//	QuadPatchFiber3d **_QuadPatchFiber3ds = 0;
//	int numQuadFibers = quadPatchFiber3ds->Length;
//	if (numQuadFibers > 0)
//	{
//		_QuadPatchFiber3ds = new QuadPatchFiber3d*[numQuadFibers];
//		for (int i = 0; i<numQuadFibers; i++)
//		{
//			_QuadPatchFiber3ds[i] = quadPatchFiber3ds[i]->_QuadPatchFiber3d;
//		}
//	}
//
//	_GJ = GJ->_UniaxialMaterial;
//	_SectionForceDeformation = new FiberSection3d(tag, numFibers, _fibers, numQuadFibers, _QuadPatchFiber3ds, _GJ);
//}
//
//FiberSection3dWrapper::FiberSection3dWrapper(int tag, FiberSectionReprWrapper ^ fiberSectionRepr,
//	Dictionary<int, UniaxialMaterialWrapper^>^ materials, array<QuadPatchFiber3dWrapper^>^ quadPatchFiber3ds, UniaxialMaterialWrapper^ GJ)
//{
//	int numPatches = fiberSectionRepr->_FiberSectionRepr->getNumPatches();
//	Patch** _patches = fiberSectionRepr->_FiberSectionRepr->getPatches();
//	int numReinfLayers = fiberSectionRepr->_FiberSectionRepr->getNumReinfLayers();
//	ReinfLayer** _reinfLayers = fiberSectionRepr->_FiberSectionRepr->getReinfLayers();
//	int numFibers = 0;
//	for (int i = 0; i < numPatches; i++)
//		numFibers += _patches[i]->getNumCells();
//	for (int i = 0; i < numReinfLayers; i++)
//		numFibers += _reinfLayers[i]->getNumReinfBars();
//
//	Fiber** _fibers = new Fiber*[numFibers];
//	int currentFiber = 0;
//
//	for (int i = 0; i < numPatches; i++)
//	{
//		int numCells = _patches[i]->getNumCells();
//		Cell** _cells = _patches[i]->getCells();
//		int theMatTag = _patches[i]->getMaterialID();
//
//		if (!materials->ContainsKey(theMatTag))
//			throw gcnew System::Exception("UniaxialMaterial with tag" + theMatTag.ToString() + " not found ?!");
//
//		UniaxialMaterial* theMat = materials[theMatTag]->_UniaxialMaterial;
//
//		for (int j = 0; j < numCells; j++)
//		{
//			_fibers[currentFiber] = new UniaxialFiber3d(currentFiber, *theMat, _cells[j]->getArea(),
//				_cells[j]->getCentroidPosition());
//			currentFiber++;
//		}
//	}
//
//
//	for (int i = 0; i < numReinfLayers; i++)
//	{
//		int numReinfBars = _reinfLayers[i]->getNumReinfBars();
//		ReinfBar* _reinfBars = _reinfLayers[i]->getReinfBars();
//		int theMatTag = _reinfLayers[i]->getMaterialID();
//
//		if (!materials->ContainsKey(theMatTag))
//			throw gcnew System::Exception("UniaxialMaterial with tag" + theMatTag.ToString() + " not found ?!");
//
//		UniaxialMaterial* theMat = materials[theMatTag]->_UniaxialMaterial;
//
//		for (int j = 0; j < numReinfBars; j++)
//		{
//			_fibers[currentFiber] = new UniaxialFiber3d(currentFiber, *theMat, _reinfBars[i].getArea(),
//				_reinfBars[i].getPosition());
//			currentFiber++;
//		}
//	}
//
//	this->_QuadPatchFiber3ds = 0;
//	numQuadFibers = quadPatchFiber3ds->Length;
//	if (numQuadFibers > 0)
//	{
//		this->_QuadPatchFiber3ds = new QuadPatchFiber3d*[numQuadFibers];
//		for (int i = 0; i<numQuadFibers; i++)
//		{
//			_QuadPatchFiber3ds[i] = quadPatchFiber3ds[i]->_QuadPatchFiber3d;
//		}
//	}
//	_GJ = GJ->_UniaxialMaterial;
//	_SectionForceDeformation = new FiberSection3d(tag, numFibers, _fibers, numQuadFibers, _QuadPatchFiber3ds, _GJ);
//}


