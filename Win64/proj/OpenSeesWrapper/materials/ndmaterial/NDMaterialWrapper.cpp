#include "stdafx.h"
#include "nDMaterialWrapper.h"

using namespace OpenSees::Materials::NDMaterials;

double* array2pointer4(array<double>^ arr) {
	double* rets = new double[arr->Length];
	for (int i = 0; i < arr->Length; i++)
	{
		rets[i] = arr[i];
	}

	return rets;
}

UniaxialMaterial** array2pointer4(array<UniaxialMaterialWrapper^>^ theMaterials) {
	UniaxialMaterial** mats = new UniaxialMaterial*[theMaterials->Length];
	for (int i = 0; i < theMaterials->Length; i++)
	{
		mats[i] = theMaterials[i]->_UniaxialMaterial;
	}

	return mats;
}

NDMaterial** array2pointer4(array<NDMaterialWrapper^>^ theMaterials) {
	NDMaterial** mats = new NDMaterial*[theMaterials->Length];
	for (int i = 0; i < theMaterials->Length; i++)
	{
		mats[i] = theMaterials[i]->_NDMaterial;
	}

	return mats;
}

NDMaterialWrapper::NDMaterialWrapper()
{

}

ElasticIsotropicMaterialWrapper::ElasticIsotropicMaterialWrapper(int tag, double E, double nu, double rho)
{
	_NDMaterial = new ElasticIsotropicMaterial(tag, E, nu, rho);
}

MultiaxialCyclicPlasticityWrapper::MultiaxialCyclicPlasticityWrapper(int    tag,
	double rho,
	double K,
	double G,
	double Su,
	double Ho_kin,
	double Parameterh,
	double Parameter_m,
	double Parameter_beta,
	double Kcoeff,
	double viscosity)
{
	_NDMaterial = new MultiaxialCyclicPlasticity(    tag,
		ND_TAG_MultiaxialCyclicPlasticity,
		 rho,
		 K,
		 G,
		 Su,
		 Ho_kin,
		 Parameterh,
		 Parameter_m,
		 Parameter_beta,
		 Kcoeff,
		 viscosity);
}

MultiaxialCyclicPlasticityWrapper::MultiaxialCyclicPlasticityWrapper(int tag, double rho, double K, double G)
{
	_NDMaterial = new MultiaxialCyclicPlasticity(tag,
		ND_TAG_MultiaxialCyclicPlasticity,
		rho,
		K,
		G);
}

MultiaxialCyclicPlasticity3DWrapper::MultiaxialCyclicPlasticity3DWrapper(int    tag,
	double rho,
	double K,
	double G,
	double Su,
	double Ho_kin,
	double Parameterh,
	double Parameter_m,
	double Parameter_beta,
	double Kcoeff,
	double viscosity)
{
	_NDMaterial = new MultiaxialCyclicPlasticity3D(tag,
		rho,
		K,
		G,
		Su,
		Ho_kin,
		Parameterh,
		Parameter_m,
		Parameter_beta,
		Kcoeff,
		viscosity);
}

MultiaxialCyclicPlasticity3DWrapper::MultiaxialCyclicPlasticity3DWrapper(int tag, double rho, double K, double G)
{
	_NDMaterial = new MultiaxialCyclicPlasticity3D( tag,  rho,  K,  G);
}

MultiaxialCyclicPlasticityAxiSymmWrapper::MultiaxialCyclicPlasticityAxiSymmWrapper(int    tag,
	double rho,
	double K,
	double G,
	double Su,
	double Ho_kin,
	double Parameterh,
	double Parameter_m,
	double Parameter_beta,
	double Kcoeff,
	double viscosity)
{
	_NDMaterial = new MultiaxialCyclicPlasticityAxiSymm(    tag,
		 rho,
		 K,
		 G,
		 Su,
		 Ho_kin,
		 Parameterh,
		 Parameter_m,
		 Parameter_beta,
		 Kcoeff,
		 viscosity);
}

MultiaxialCyclicPlasticityAxiSymmWrapper::MultiaxialCyclicPlasticityAxiSymmWrapper(int tag, double rho, double K, double G)
{
	_NDMaterial = new MultiaxialCyclicPlasticityAxiSymm(tag, rho, K, G);
}


MultiaxialCyclicPlasticityPlaneStrainWrapper::MultiaxialCyclicPlasticityPlaneStrainWrapper(int    tag,
	double rho,
	double K,
	double G,
	double Su,
	double Ho_kin,
	double Parameterh,
	double Parameter_m,
	double Parameter_beta,
	double Kcoeff,
	double viscosity)
{
	_NDMaterial = new MultiaxialCyclicPlasticityPlaneStrain(tag,
		rho,
		K,
		G,
		Su,
		Ho_kin,
		Parameterh,
		Parameter_m,
		Parameter_beta,
		Kcoeff,
		viscosity);
}

MultiaxialCyclicPlasticityPlaneStrainWrapper::MultiaxialCyclicPlasticityPlaneStrainWrapper(int tag, double rho, double K, double G)
{
	_NDMaterial = new MultiaxialCyclicPlasticityPlaneStrain(tag, rho, K, G);
}

ElasticIsotropicMaterialThermalWrapper::ElasticIsotropicMaterialThermalWrapper(int tag, double E, double nu, double rho, double alpha, int softindex)
{
	_NDMaterial = new ElasticIsotropicMaterialThermal(tag, E, nu, rho, alpha, softindex);
}

ElasticIsotropic3DThermalWrapper::ElasticIsotropic3DThermalWrapper(int tag, double E, double nu, double rho, double alpha, int softindex)
{
	_NDMaterial = new ElasticIsotropic3DThermal(tag, E, nu, rho, alpha, softindex);
}

ElasticIsotropicAxiSymmWrapper::ElasticIsotropicAxiSymmWrapper(int tag, double E, double nu, double rho)
{
	_NDMaterial = new ElasticIsotropicAxiSymm(tag, E, nu, rho);
}

ElasticIsotropicBeamFiberWrapper::ElasticIsotropicBeamFiberWrapper(int tag, double E, double nu, double rho)
{
	_NDMaterial = new ElasticIsotropicBeamFiber(tag, E, nu, rho);
}

ElasticIsotropicBeamFiber2dWrapper::ElasticIsotropicBeamFiber2dWrapper(int tag, double E, double nu, double rho)
{
	_NDMaterial = new ElasticIsotropicBeamFiber2d(tag, E, nu, rho);
}

ElasticIsotropicPlaneStress2DWrapper::ElasticIsotropicPlaneStress2DWrapper(int tag, double E, double nu, double rho)
{
	_NDMaterial = new ElasticIsotropicPlaneStress2D(tag, E, nu, rho);
}

ElasticIsotropicPlaneStrain2DWrapper::ElasticIsotropicPlaneStrain2DWrapper(int tag, double E, double nu, double rho)
{
	_NDMaterial = new ElasticIsotropicPlaneStrain2D(tag, E, nu, rho);
}

ElasticIsotropicPlateFiberWrapper::ElasticIsotropicPlateFiberWrapper(int tag, double E, double nu, double rho)
{
	_NDMaterial = new ElasticIsotropicPlateFiber(tag, E, nu, rho);
}

ElasticIsotropicThreeDimensionalWrapper::ElasticIsotropicThreeDimensionalWrapper(int tag, double E, double nu, double rho)
{
	_NDMaterial = new ElasticIsotropicThreeDimensional(tag, E, nu, rho);
}

PressureDependentElastic3DWrapper::PressureDependentElastic3DWrapper(int tag,
	double E,
	double nu,
	double rhop,
	double expp,
	double pr,
	double pop)
{
	_NDMaterial = new PressureDependentElastic3D(tag,
		E,
		nu,
		rhop,
		expp,
		pr,
		pop);
}

FeapMaterialWrapper::FeapMaterialWrapper(int tag, int classTag, int numHV, int numData,
	double rho)
{
	_NDMaterial = new FeapMaterial(tag, classTag, numHV, numData,
		rho);
}

FeapMaterial01Wrapper::FeapMaterial01Wrapper(int tag, double E, double nu, double rho)
{
	_NDMaterial = new FeapMaterial01( tag,  E,  nu,  rho);
}

FeapMaterial02Wrapper::FeapMaterial02Wrapper(int tag, double K, double G, double muK, double muG,
	double lamK, double lamG, double theta)
{
	_NDMaterial = new FeapMaterial02( tag,  K,  G,  muK,  muG,
		 lamK,  lamG,  theta);
}

FeapMaterial03Wrapper::FeapMaterial03Wrapper(int tag, double K, double G, double sigY, double Hiso)
{
	_NDMaterial = new FeapMaterial03(tag, K, G, sigY, Hiso);
}

J2PlasticityWrapper::J2PlasticityWrapper(int    tag,
	int    classTag,
	double K,
	double G,
	double yield0,
	double yield_infty,
	double d,
	double H,
	double viscosity,
	double rho)
{
	_NDMaterial = new J2Plasticity(tag,
		classTag,
		K,
		G,
		yield0,
		yield_infty,
		d,
		H,
		viscosity,
		rho);
}

J2PlasticityWrapper::J2PlasticityWrapper(int    tag,
	int    classTag,
	double K,
	double G
	)
{
	_NDMaterial = new J2Plasticity(tag,
		classTag,
		K,
		G
	);
}

J2AxiSymmWrapper::J2AxiSymmWrapper(int    tag,
	double K,
	double G,
	double yield0,
	double yield_infty,
	double d,
	double H,
	double viscosity,
	double rho
)
{
	_NDMaterial = new J2AxiSymm(tag,
		K,
		G,
		yield0,
		yield_infty,
		d,
		H,
		viscosity,
		rho
	);
}

J2AxiSymmWrapper::J2AxiSymmWrapper(int tag, double K, double G
)
{
	_NDMaterial = new J2AxiSymm(tag,
		K,
		G);
}

J2BeamFiber2dWrapper::J2BeamFiber2dWrapper(int tag, double E, double G, double sigY, double Hi, double Hk)
{
	_NDMaterial = new J2BeamFiber2d(tag, E, G, sigY, Hi, Hk);
}

J2BeamFiber3dWrapper::J2BeamFiber3dWrapper(int tag, double E, double G, double sigY, double Hi, double Hk)
{
	_NDMaterial = new J2BeamFiber3d(tag, E, G, sigY, Hi, Hk);
}

J2CyclicBoundingSurfaceWrapper::J2CyclicBoundingSurfaceWrapper(int    tag,
	double G,
	double K,
	double su,
	double rho,
	double h,
	double m,
	double k_in,
	double beta)
{
	_NDMaterial = new J2CyclicBoundingSurface(tag,
		G,
		K,
		su,
		rho,
		h,
		m,
		k_in,
		beta);
}

J2PlaneStrainWrapper::J2PlaneStrainWrapper(int    tag,
	double K,
	double G,
	double yield0,
	double yield_infty,
	double d,
	double H,
	double viscosity)
{
	_NDMaterial = new J2PlaneStrain(tag,
		K,
		G,
		yield0,
		yield_infty,
		d,
		H,
		viscosity);
}

J2PlaneStressWrapper::J2PlaneStressWrapper(int    tag,
	double K,
	double G,
	double yield0,
	double yield_infty,
	double d,
	double H,
	double viscosity,
	double rho)
{
	_NDMaterial = new J2PlaneStress(tag,
		K,
		G,
		yield0,
		yield_infty,
		d,
		H,
		viscosity,
		rho);
}

J2PlaneStressWrapper::J2PlaneStressWrapper(int tag, double K, double G)
{
	_NDMaterial = new J2PlaneStress(tag, K, G);
}

J2PlasticityThermalWrapper::J2PlasticityThermalWrapper(int    tag,
	int    classTag,
	double K,
	double G,
	double yield0,
	double yield_infty,
	double d,
	double H,
	double viscosity,
	double rho) {
	_NDMaterial = new J2PlasticityThermal(tag, classTag, K, G, yield0, yield_infty, d, H, viscosity, rho);
}

J2PlasticityThermalWrapper::J2PlasticityThermalWrapper(int tag, int classTag, double K, double G) {
	_NDMaterial = new J2PlasticityThermal(tag, classTag, K, G);
}

J2PlateFiberWrapper::J2PlateFiberWrapper(int    tag,
	double K,
	double G,
	double yield0,
	double yield_infty,
	double d,
	double H,
	double viscosity,
	double rho) {
	_NDMaterial = new J2PlateFiber(tag, K, G, yield0, yield_infty, d, H, viscosity, rho);
}

J2PlateFiberWrapper::J2PlateFiberWrapper(int tag, double K, double G) {
	_NDMaterial = new J2PlateFiber(tag, K, G);
}

J2PlateFibreWrapper::J2PlateFibreWrapper(int tag, double E, double G, double sigY, double Hi, double Hk) {
	_NDMaterial = new J2PlateFibre(tag, E, G, sigY, Hi, Hk);
}

J2ThreeDimensionalWrapper::J2ThreeDimensionalWrapper(int    tag,
	double K,
	double G,
	double yield0,
	double yield_infty,
	double d,
	double H,
	double viscosity,
	double rho) {
	_NDMaterial = new J2ThreeDimensional(tag, K, G, yield0, yield_infty, d, H, viscosity, rho);
}

J2ThreeDimensionalThermalWrapper::J2ThreeDimensionalThermalWrapper(int    tag,
	double K,
	double G,
	double yield0,
	double yield_infty,
	double d,
	double H,
	double viscosity,
	double rho) {
	_NDMaterial = new J2ThreeDimensionalThermal(tag, K, G, yield0, yield_infty, d, H, viscosity, rho);
}

J2ThreeDimensionalThermalWrapper::J2ThreeDimensionalThermalWrapper(int tag, double K, double G) {
	_NDMaterial = new J2ThreeDimensionalThermal(tag, K, G);
}

MaterialCMMWrapper::MaterialCMMWrapper(int tag, int layer, array<double>^ props) {

	_NDMaterial = new MaterialCMM(tag, layer, array2pointer4(props));
}


FAFourSteelPCPlaneStressWrapper::FAFourSteelPCPlaneStressWrapper(int      tag,
	double   RHO,
	UniaxialMaterialWrapper^ t1,
	UniaxialMaterialWrapper^ t2,
	UniaxialMaterialWrapper^ s1,
	UniaxialMaterialWrapper^ s2,
	UniaxialMaterialWrapper^ c1,
	UniaxialMaterialWrapper^ c2,
	double   ANGLE1,
	double   ANGLE2,
	double   ANGLE3,
	double   ANGLE4,
	double   ROU1,
	double   ROU2,
	double   ROU3,
	double   ROU4,
	double	PSTRAIN1,
	double	PSTRAIN2,
	double   FPC,
	double   FY1,
	double	FY2,
	double   E,
	double   EPSC0) {
	_NDMaterial = new FAFourSteelPCPlaneStress(tag, RHO, t1->_UniaxialMaterial, t2->_UniaxialMaterial, s1->_UniaxialMaterial, s2->_UniaxialMaterial, c1->_UniaxialMaterial, c2->_UniaxialMaterial, ANGLE1, ANGLE2, ANGLE3, ANGLE4, ROU1, ROU2, ROU3, ROU4,
			PSTRAIN1,
			PSTRAIN2, FPC, FY1,
			FY2, E, EPSC0);
}

FAFourSteelRCPlaneStressWrapper::FAFourSteelRCPlaneStressWrapper(int      tag,
	double   RHO,
	UniaxialMaterialWrapper^ t1,
	UniaxialMaterialWrapper^ t2,
	UniaxialMaterialWrapper^ s1,
	UniaxialMaterialWrapper^ s2,
	UniaxialMaterialWrapper^ c1,
	UniaxialMaterialWrapper^ c2,
	double   ANGLE1,
	double   ANGLE2,
	double   ANGLE3,
	double   ANGLE4,
	double   ROU1,
	double   ROU2,
	double   ROU3,
	double   ROU4,
	double   FPC,
	double   FY,
	double   E,
	double   EPSC0) {
	_NDMaterial = new FAFourSteelRCPlaneStress(tag, RHO, t1->_UniaxialMaterial, t2->_UniaxialMaterial, s1->_UniaxialMaterial, s2->_UniaxialMaterial, c1->_UniaxialMaterial, c2->_UniaxialMaterial, ANGLE1, ANGLE2, ANGLE3, ANGLE4, ROU1, ROU2, ROU3, ROU4, FPC, FY, E, EPSC0);
}

FAPrestressedConcretePlaneStressWrapper::FAPrestressedConcretePlaneStressWrapper(int      tag,
	double   RHO,
	UniaxialMaterialWrapper^ t1,
	UniaxialMaterialWrapper^ s1,
	UniaxialMaterialWrapper^ c1,
	UniaxialMaterialWrapper^ c2,
	double   ANGLE1,
	double   ANGLE2,
	double   ROU1,
	double   ROU2,
	double	PSTRAIN,
	double   FPC,
	double   FY1,
	double	FY2,
	double   E,
	double   EPSC0) {
	_NDMaterial = new FAPrestressedConcretePlaneStress(tag, RHO, t1->_UniaxialMaterial, s1->_UniaxialMaterial, c1->_UniaxialMaterial, c2->_UniaxialMaterial, ANGLE1, ANGLE2, ROU1, ROU2,
			PSTRAIN, FPC, FY1,
			FY2, E, EPSC0);
}


FAReinforcedConcretePlaneStressWrapper::FAReinforcedConcretePlaneStressWrapper(int      tag,
	double   RHO,
	UniaxialMaterialWrapper^ t1,
	UniaxialMaterialWrapper^ s1,
	UniaxialMaterialWrapper^ c1,
	UniaxialMaterialWrapper^ c2,
	double   ANGLE1,
	double   ANGLE2,
	double   ROU1,
	double   ROU2,
	double   FPC,
	double   FY,
	double   E,
	double   EPSC0) {
	_NDMaterial = new FAReinforcedConcretePlaneStress(tag, RHO, t1->_UniaxialMaterial, s1->_UniaxialMaterial, c1->_UniaxialMaterial, c2->_UniaxialMaterial, ANGLE1, ANGLE2, ROU1, ROU2, FPC, FY, E, EPSC0);
}


PrestressedConcretePlaneStressWrapper::PrestressedConcretePlaneStressWrapper(int      tag,
	double   RHO,
	UniaxialMaterialWrapper^ t1,
	UniaxialMaterialWrapper^ s1,
	UniaxialMaterialWrapper^ c1,
	UniaxialMaterialWrapper^ c2,
	double   ANGLE1,
	double   ANGLE2,
	double   ROU1,
	double   ROU2,
	double	PSTRAIN,
	double   FPC,
	double   FY1,
	double	FY2,
	double   E,
	double   EPSC0) {
	_NDMaterial = new PrestressedConcretePlaneStress(tag, RHO, t1->_UniaxialMaterial, s1->_UniaxialMaterial, c1->_UniaxialMaterial, c2->_UniaxialMaterial, ANGLE1, ANGLE2, ROU1, ROU2, PSTRAIN, FPC, FY1, FY2, E, EPSC0);
}

RAFourSteelPCPlaneStressWrapper::RAFourSteelPCPlaneStressWrapper(int      tag,
	double   RHO,
	UniaxialMaterialWrapper^ t1,
	UniaxialMaterialWrapper^ t2,
	UniaxialMaterialWrapper^ s1,
	UniaxialMaterialWrapper^ s2,
	UniaxialMaterialWrapper^ c1,
	UniaxialMaterialWrapper^ c2,
	double   ANGLE1,
	double   ANGLE2,
	double   ANGLE3,
	double   ANGLE4,
	double   ROU1,
	double   ROU2,
	double   ROU3,
	double   ROU4,
	double	PSTRAIN1,
	double	PSTRAIN2,
	double   FPC,
	double   FY1,
	double	FY2,
	double   E,
	double   EPSC0) {
	_NDMaterial = new RAFourSteelPCPlaneStress(tag, RHO, t1->_UniaxialMaterial, t2->_UniaxialMaterial, s1->_UniaxialMaterial, s2->_UniaxialMaterial, c1->_UniaxialMaterial, c2->_UniaxialMaterial, ANGLE1, ANGLE2, ANGLE3, ANGLE4, ROU1, ROU2, ROU3, ROU4, PSTRAIN1, PSTRAIN2, FPC, FY1, FY2, E, EPSC0);
}

RAFourSteelRCPlaneStressWrapper::RAFourSteelRCPlaneStressWrapper(int      tag,
	double   RHO,
	UniaxialMaterialWrapper^ s1,
	UniaxialMaterialWrapper^ s2,
	UniaxialMaterialWrapper^ s3,
	UniaxialMaterialWrapper^ s4,
	UniaxialMaterialWrapper^ c1,
	UniaxialMaterialWrapper^ c2,
	double   ANGLE1,
	double   ANGLE2,
	double   ANGLE3,
	double   ANGLE4,
	double   ROU1,
	double   ROU2,
	double   ROU3,
	double   ROU4,
	double   FPC,
	double   FY,
	double   E,
	double   EPSC0) {
	_NDMaterial = new RAFourSteelRCPlaneStress(tag, RHO, s1->_UniaxialMaterial, s2->_UniaxialMaterial, s3->_UniaxialMaterial, s4->_UniaxialMaterial, c1->_UniaxialMaterial, c2->_UniaxialMaterial, ANGLE1, ANGLE2, ANGLE3, ANGLE4, ROU1, ROU2, ROU3, ROU4, FPC, FY, E, EPSC0);
}


ReinforcedConcretePlaneStressWrapper::ReinforcedConcretePlaneStressWrapper(int      tag,
	double   RHO,
	UniaxialMaterialWrapper^ s1,
	UniaxialMaterialWrapper^ s2,
	UniaxialMaterialWrapper^ c1,
	UniaxialMaterialWrapper^ c2,
	double   ANGLE1,
	double   ANGLE2,
	double   ROU1,
	double   ROU2,
	double   FPC,
	double   FY,
	double   E,
	double   EPSC0) {
	_NDMaterial = new ReinforcedConcretePlaneStress(tag, RHO, s1->_UniaxialMaterial, s2->_UniaxialMaterial, c1->_UniaxialMaterial, c2->_UniaxialMaterial, ANGLE1, ANGLE2, ROU1, ROU2, FPC, FY, E, EPSC0);
}

FluidSolidPorousMaterialWrapper::FluidSolidPorousMaterialWrapper(int tag, int nd, NDMaterialWrapper^ soilMat,
	double combinedBulkModul, double atm) {
	_NDMaterial = new FluidSolidPorousMaterial(tag, nd, *soilMat->_NDMaterial, combinedBulkModul, atm);
}

MultiYieldSurfaceClayWrapper::MultiYieldSurfaceClayWrapper(int tag,
	int nd,
	double rho,
	double refShearModul,
	double refBulkModul,
	double cohesi,
	double peakShearStra,
	double frictionAng,
	double refPress,
	double pressDependCoe,
	int   numberOfYieldSurf,
	array<double>^ gredu) {
	_NDMaterial = new MultiYieldSurfaceClay(tag, nd, rho, refShearModul, refBulkModul, cohesi, peakShearStra, frictionAng, refPress, pressDependCoe,
		numberOfYieldSurf, array2pointer4(gredu));
}


PressureDependMultiYieldWrapper::PressureDependMultiYieldWrapper(int tag,
	int nd,
	double rho,
	double refShearModul,
	double refBulkModul,
	double frictionAng,
	double peakShearStra,
	double refPress,
	double pressDependCoe,
	double phaseTransformAngle,
	double contractionParam1,
	double dilationParam1,
	double dilationParam2,
	double liquefactionParam1,
	double liquefactionParam2,
	double liquefactionParam4,
	int   numberOfYieldSurf,
	array<double>^ gredu,
	double e,
	double volLimit1,
	double volLimit2,
	double volLimit3,
	double atm,
	double cohesi,
	double hv,
	double pv) {
	_NDMaterial = new PressureDependMultiYield(tag, nd, rho, refShearModul, refBulkModul, frictionAng, peakShearStra,
		refPress, pressDependCoe, phaseTransformAngle, contractionParam1, dilationParam1, dilationParam2,
		liquefactionParam1, liquefactionParam2, liquefactionParam4, numberOfYieldSurf, array2pointer4(gredu), e, volLimit1, volLimit2, volLimit3, atm, cohesi, hv, pv);
}

PressureDependMultiYield02Wrapper::PressureDependMultiYield02Wrapper(int tag,
	int nd,
	double rho,
	double refShearModul,
	double refBulkModul,
	double frictionAng,
	double peakShearStra,
	double refPress,
	double pressDependCoe,
	double phaseTransformAngle,
	double contractionParam1,
	double contractionParam3,
	double dilationParam1,
	double dilationParam3,
	int   numberOfYieldSurf,
	array<double>^ gredu,
	double contractionParam2,
	double dilationParam2,
	double liquefactionParam1,
	double liquefactionParam2,
	double e,
	double volLimit1,
	double volLimit2,
	double volLimit3,
	double atm,
	double cohesi,
	double hv,
	double pv) {
	_NDMaterial = new PressureDependMultiYield02(tag, nd, rho, refShearModul, refBulkModul, frictionAng, peakShearStra, refPress, pressDependCoe,
		phaseTransformAngle, contractionParam1, contractionParam3, dilationParam1, dilationParam3, numberOfYieldSurf, array2pointer4(gredu),
		contractionParam2, dilationParam2, liquefactionParam1, liquefactionParam2, e, volLimit1, volLimit2, volLimit3, atm, cohesi, hv, pv);
}

PressureIndependMultiYieldWrapper::PressureIndependMultiYieldWrapper(int tag,
	int nd,
	double rho,
	double refShearModul,
	double refBulkModul,
	double cohesi,
	double peakShearStra,
	double frictionAng,
	double refPress,
	double pressDependCoe,
	int   numberOfYieldSurf,
	array<double>^ gredu) {
	_NDMaterial = new PressureIndependMultiYield(tag, nd, rho, refShearModul, refBulkModul, cohesi, peakShearStra, frictionAng, refPress,
		pressDependCoe, numberOfYieldSurf, array2pointer4(gredu));
}

//StressDensityModelWrapper::StressDensityModelWrapper(int tag, int classTag, double constDensity,
//	// SD model  parameters		
//	double initialVoidRatio, double constA, double exponentN,
//	double poissonRatio, double constAlpha1, double constBeta1,
//	double constAlpha2, double constBeta2, double constAlpha3,
//	double constBeta3, double constDegradation, double constMumin,
//	double constMucyclic, double constDilatancyStrain,
//	double constMumax, double constPatm) {
//	_NDMaterial = new StressDensityModel(tag, classTag, constDensity, initialVoidRatio, constA, exponentN, poissonRatio, constAlpha1, constBeta1, constAlpha2, constBeta2, constAlpha3, constBeta3, constDegradation, constMumin, constMucyclic, constDilatancyStrain, constMumax, constPatm);
//}
//
//
//StressDensityModelWrapper::StressDensityModelWrapper(int tag, int classTag, double constDensity,
//	// SD model  parameters		
//	double initialVoidRatio, double constA, double exponentN,
//	double poissonRatio, double constAlpha1, double constBeta1,
//	double constAlpha2, double constBeta2, double constAlpha3,
//	double constBeta3, double constDegradation, double constMumin,
//	double constMucyclic, double constDilatancyStrain,
//	double constMumax, double constPatm, double constsslvoidatP1, double constsslvoidatP2, double constsslvoidatP3,
//	double constsslvoidatP4, double constsslvoidatP5, double constsslvoidatP6,
//	double constsslvoidatP7, double constsslvoidatP8, double constsslvoidatP9,
//	double constsslvoidatP10,
//	// hydrostatic state line void ratio
//	double consthslvoid,
//	// reference pressures
//	double constP1, double constP2, double constP3, double constP4,
//	double constP5, double constP6, double constP7, double constP8,
//	double constP9, double constP10,
//	// offset of the failure surface
//	double constRxx, double constRyy, double constRzz,
//	double constRxy, double constRyz, double constRzx) {
//	_NDMaterial = new StressDensityModel(tag, classTag, constDensity, initialVoidRatio, constA, exponentN, poissonRatio, constAlpha1, constBeta1, constAlpha2, constBeta2, constAlpha3, constBeta3, constDegradation, constMumin, constMucyclic, constDilatancyStrain, constMumax, constPatm, constsslvoidatP1, constsslvoidatP2, constsslvoidatP3, constsslvoidatP4, constsslvoidatP5, constsslvoidatP6, constsslvoidatP7, constsslvoidatP8, constsslvoidatP9, constsslvoidatP10, consthslvoid, constP1, constP2, constP3, constP4, constP5, constP6, constP7, constP8, constP9, constP10, constRxx, constRyy, constRzz, constRxy, constRyz, constRzx);
//}
//
//StressDensityModel2DWrapper::StressDensityModel2DWrapper(int tag, double constDensity,
//	// SD model  parameters
//	double initialVoidRatio, double constA, double exponentN, double poissonRatio,
//	double constAlpha1, double constBeta1, double constAlpha2, double constBeta2,
//	double constAlpha3, double constBeta3, double constDegradation, double constMumin,
//	double constMucyclic, double constDilatancyStrain, double constMumax, double constPatm,
//	// steady state line void ratio
//	double constsslvoidatP1, double constsslvoidatP2, double constsslvoidatP3,
//	double constsslvoidatP4, double constsslvoidatP5, double constsslvoidatP6,
//	double constsslvoidatP7, double constsslvoidatP8, double constsslvoidatP9,
//	double constsslvoidatP10,
//	// hydrostatic state line void ratio
//	double consthslvoid,
//	// reference pressures 
//	double constP1, double constP2, double constP3, double constP4, double constP5,
//	double constP6, double constP7, double constP8, double constP9, double constP10,
//	// offset of the failure surface
//	double constRxx, double constRyy, double constRxy) {
//	_NDMaterial = new StressDensityModel2D(tag, constDensity, initialVoidRatio, constA, exponentN, poissonRatio, constAlpha1, constBeta1, constAlpha2, constBeta2, constAlpha3, constBeta3, constDegradation, constMumin, constMucyclic, constDilatancyStrain, constMumax, constPatm, constsslvoidatP1, constsslvoidatP2, constsslvoidatP3, constsslvoidatP4, constsslvoidatP5, constsslvoidatP6, constsslvoidatP7, constsslvoidatP8, constsslvoidatP9, constsslvoidatP10, consthslvoid, constP1, constP2, constP3, constP4, constP5, constP6, constP7, constP8, constP9, constP10, constRxx, constRyy, constRxy);
//}
//
//StressDensityModel3DWrapper::StressDensityModel3DWrapper(int tag, double constDensity,
//	// SD model  parameters
//	double initialVoidRatio, double constA, double exponentN, double poissonRatio,
//	double constAlpha1, double constBeta1, double constAlpha2, double constBeta2,
//	double constAlpha3, double constBeta3, double constDegradation, double constMumin,
//	double constMucyclic, double constDilatancyStrain, double constMumax, double constPatm,
//	// steady state line void ratio
//	double constsslvoidatP1, double constsslvoidatP2, double constsslvoidatP3,
//	double constsslvoidatP4, double constsslvoidatP5, double constsslvoidatP6,
//	double constsslvoidatP7, double constsslvoidatP8, double constsslvoidatP9,
//	double constsslvoidatP10,
//	// hydrostatic state line void ratio
//	double consthslvoid,
//	// reference pressures 
//	double constP1, double constP2, double constP3, double constP4, double constP5,
//	double constP6, double constP7, double constP8, double constP9, double constP10,
//	// offset of the failure surface
//	double constRxx, double constRyy, double constRzz, double constRxy,
//	double constRyz, double constRzx) {
//	_NDMaterial = new StressDensityModel3D(tag, constDensity, initialVoidRatio, constA, exponentN, poissonRatio, constAlpha1, constBeta1, constAlpha2, constBeta2, constAlpha3, constBeta3, constDegradation, constMumin, constMucyclic, constDilatancyStrain, constMumax, constPatm, constsslvoidatP1, constsslvoidatP2, constsslvoidatP3, constsslvoidatP4, constsslvoidatP5, constsslvoidatP6, constsslvoidatP7, constsslvoidatP8, constsslvoidatP9, constsslvoidatP10, consthslvoid, constP1, constP2, constP3, constP4, constP5, constP6, constP7, constP8, constP9, constP10, constRxx, constRyy, constRzz, constRxy, constRyz, constRzx);
//}

BeamFiberMaterialWrapper::BeamFiberMaterialWrapper(int tag, NDMaterialWrapper^ theMat) {
	_NDMaterial = new BeamFiberMaterial(tag, *theMat->_NDMaterial);
}


BeamFiberMaterial2dWrapper::BeamFiberMaterial2dWrapper(int tag, NDMaterialWrapper^ theMat) {
	_NDMaterial = new BeamFiberMaterial2d(tag, *theMat->_NDMaterial);
}

PlaneStrainMaterialWrapper::PlaneStrainMaterialWrapper(int tag, NDMaterialWrapper^ theMat) {
	_NDMaterial = new PlaneStrainMaterial(tag, *theMat->_NDMaterial);
}

PlaneStressLayeredMaterialWrapper::PlaneStressLayeredMaterialWrapper(int tag,
	int iLayers,
	array<double>^ thickness,
	array<NDMaterialWrapper^>^ fibers) {
	_NDMaterial = new PlaneStressLayeredMaterial(tag, iLayers, array2pointer4(thickness), array2pointer4(fibers));
}

PlaneStressMaterialWrapper::PlaneStressMaterialWrapper(int tag, NDMaterialWrapper^ the3DMaterial) {
	_NDMaterial = new PlaneStressMaterial(tag, *the3DMaterial->_NDMaterial);
}

PlaneStressRebarMaterialWrapper::PlaneStressRebarMaterialWrapper(int tag, UniaxialMaterialWrapper^ uniMat,
	double ang) {
	_NDMaterial = new PlaneStressRebarMaterial(tag, *uniMat->_UniaxialMaterial, ang);
}

PlaneStressSimplifiedJ2Wrapper::PlaneStressSimplifiedJ2Wrapper(int tag, int nd, NDMaterialWrapper^ uniMat) {
	_NDMaterial = new PlaneStressSimplifiedJ2(tag, nd, *uniMat->_NDMaterial);
}

PlaneStressUserMaterialWrapper::PlaneStressUserMaterialWrapper(int tag, int istatevs, int iprops, array<double>^ props) {
	_NDMaterial = new PlaneStressUserMaterial(tag, istatevs, iprops, array2pointer4(props));
}

PlateFiberMaterialWrapper::PlateFiberMaterialWrapper(int    tag,
	NDMaterialWrapper^ the3DMaterial) {
	_NDMaterial = new PlateFiberMaterial(tag, *the3DMaterial->_NDMaterial);
}

PlateFiberMaterialThermalWrapper::PlateFiberMaterialThermalWrapper(int    tag,
	NDMaterialWrapper^ the3DMaterial) {
	_NDMaterial = new PlateFiberMaterialThermal(tag, *the3DMaterial->_NDMaterial);
}

PlateFromPlaneStressMaterialWrapper::PlateFromPlaneStressMaterialWrapper(int    tag,
	NDMaterialWrapper^ ndMat, double g) {
	_NDMaterial = new PlateFromPlaneStressMaterial(tag, *ndMat->_NDMaterial, g);
}

PlateFromPlaneStressMaterialThermalWrapper::PlateFromPlaneStressMaterialThermalWrapper(int    tag,
	NDMaterialWrapper^ ndMat, double g) {
	_NDMaterial = new PlateFromPlaneStressMaterialThermal(tag, *ndMat->_NDMaterial, g);
}

PlateRebarMaterialWrapper::PlateRebarMaterialWrapper(int    tag,
	UniaxialMaterialWrapper^ uniMat, double ang) {
	_NDMaterial = new PlateRebarMaterial(tag, *uniMat->_UniaxialMaterial, ang);
}

PlateRebarMaterialThermalWrapper::PlateRebarMaterialThermalWrapper(int    tag,
	UniaxialMaterialWrapper^ uniMat, double ang) {
	_NDMaterial = new PlateRebarMaterialThermal(tag, *uniMat->_UniaxialMaterial, ang);
}

CycLiqCPSPWrapper::CycLiqCPSPWrapper(int    tag,
	int classTag,
	double G01,
	double kappa1,
	double h1,
	double Mfc1,       //critical state
	double dre11,
	double dre21,
	double rdr1,
	double eta1,
	double dir1,
	double lamdac1,
	double ksi1,
	double e01,
	double nb1,
	double nd1,
	double ein1,      //initial void ratio
	double rho1) {
	_NDMaterial = new CycLiqCPSP(tag, classTag, G01, kappa1, h1, Mfc1, dre11, dre21, rdr1, eta1, dir1, lamdac1, ksi1, e01, nb1, nd1, ein1, rho1);
}

CycLiqCPWrapper::CycLiqCPWrapper(int    tag,
	int classTag,
	double G01,
	double kappa1,
	double h1,
	double Mfc1,       //critical state
	double dre11,
	double Mdc1,
	double dre21,
	double rdr1,
	double eta1,
	double dir1,
	double ein1,      //initial void ratio
	double rho1) {
	_NDMaterial = new CycLiqCP(tag, classTag, G01, kappa1, h1, Mfc1, dre11, Mdc1, dre21, rdr1, eta1, dir1, ein1, rho1);
}

CycLiqCPSP3DWrapper::CycLiqCPSP3DWrapper(int    tag,
	int classTag,
	double G01,
	double kappa1,
	double h1,
	double Mfc1,       //critical state
	double dre11,
	double dre21,
	double rdr1,
	double eta1,
	double dir1,
	double lamdac1,
	double ksi1,
	double e01,
	double nb1,
	double nd1,
	double ein1,      //initial void ratio
	double rho1)  {
	_NDMaterial = new CycLiqCPSP3D(tag, classTag, G01, kappa1, h1, Mfc1, dre11, dre21, rdr1, eta1, dir1, lamdac1, ksi1, e01, nb1, nd1, ein1, rho1);
}

CycLiqCPPlaneStrainWrapper::CycLiqCPPlaneStrainWrapper(int    tag,
	double G01,
	double kappa1,
	double h1,
	double Mfc1,       //critical state
	double dre11,
	double Mdc1,
	double dre21,
	double rdr1,
	double eta1,
	double dir1,
	double ein1,      //initial void ratio
	double rho1) {
	_NDMaterial = new CycLiqCPPlaneStrain(tag, G01, kappa1, h1, Mfc1, dre11, Mdc1, dre21, rdr1, eta1, dir1, ein1, rho1);
}

CycLiqCP3DWrapper::CycLiqCP3DWrapper(int    tag,
	double G01,
	double kappa1,
	double h1,
	double Mfc1,       //critical state
	double dre11,
	double Mdc1,
	double dre21,
	double rdr1,
	double eta1,
	double dir1,
	double ein1,      //initial void ratio
	double rho1) {
	_NDMaterial = new CycLiqCP3D(tag, G01, kappa1, h1, Mfc1, dre11, Mdc1, dre21, rdr1, eta1, dir1, ein1, rho1);
}

ContactMaterial3DWrapper::ContactMaterial3DWrapper(int tag, double mu, double G, double c, double t) {
	_NDMaterial = new ContactMaterial3D(tag, mu, G, c, t);
}

ContactMaterial2DWrapper::ContactMaterial2DWrapper(int tag, double mu, double G, double c, double t) {
	_NDMaterial = new ContactMaterial2D(tag, mu, G, c, t);
}

ConcreteSWrapper::ConcreteSWrapper(int tag, double rE, double rnu, double rfc, double rft, double rEs) {
	_NDMaterial = new ConcreteS(tag, rE, rnu, rfc, rft, rEs);
}

CapPlasticityWrapper::CapPlasticityWrapper(int    tag,
	double G,
	double K,
	double rho,
	double X,
	double D,
	double W,
	double R,
	double lambda,
	double theta,
	double beta,
	double alpha,
	double T,
	int ndm,
	double pTol_k) {
	_NDMaterial = new CapPlasticity(tag, G, K, rho, X, D, W, R, lambda, theta, beta, alpha, T, ndm, pTol_k);
}

CycLiqCPSPPlaneStrainWrapper::CycLiqCPSPPlaneStrainWrapper(int    tag,
	double G01,
	double kappa1,
	double h1,
	double Mfc1,       //critical state
	double dre11,
	double Mdc1,
	double dre21,
	double rdr1,
	double eta1,
	double dir1,
	double lamdac1,
	double ksi1,
	double e01,
	double nb1,
	double nd1,
	double ein1,      //initial void ratio
	double rho1) {
	_NDMaterial = new CycLiqCPSPPlaneStrain(tag, G01, kappa1, h1, Mfc1, dre11, Mdc1, dre21, rdr1, eta1, dir1, lamdac1, ksi1, e01, nb1, nd1, ein1, rho1);
}

DruckerPragerWrapper::DruckerPragerWrapper(int tag, int classTag, double bulk, double shear,
	double s_y, double r, double r_bar, double Kinfinity, double Kinit,
	double d1, double d2, double H, double t, double massDen, double atm) {
	_NDMaterial = new DruckerPrager(tag, classTag, bulk, shear, s_y, r, r_bar, Kinfinity, Kinit, d1, d2, H, t, massDen, atm);
}

DruckerPrager3DThermalWrapper::DruckerPrager3DThermalWrapper(int tag, double bulk, double shear,
	double s_y, double r, double r_bar, double Kinfinity, double Kinit,
	double d1, double d2, double H, double t, double massDen, double atm) {
	_NDMaterial = new DruckerPrager3DThermal(tag, bulk, shear, s_y, r, r_bar, Kinfinity, Kinit, d1, d2, H, t, massDen, atm);
}

DruckerPragerPlaneStrainWrapper::DruckerPragerPlaneStrainWrapper(int tag, double bulk, double shear,
	double s_y, double r, double r_bar, double Kinfinity, double Kinit,
	double d1, double d2, double H, double t, double massDens, double atm) {
	_NDMaterial = new DruckerPragerPlaneStrain(tag, bulk, shear, s_y, r, r_bar, Kinfinity, Kinit, d1, d2, H, t, massDens, atm);
}

ElasticOrthotropicMaterialWrapper::ElasticOrthotropicMaterialWrapper(int tag, int classTag,
	double Ex, double Ey, double Ez,
	double vxy, double vyz, double vzx,
	double Gxy, double Gyz, double Gz, double rho) {
	_NDMaterial = new ElasticOrthotropicMaterial(tag, classTag, Ex, Ey, Ez, vxy, vyz, vzx, Gxy, Gyz, Gz, rho);
}

ElasticOrthotropicThreeDimensionalWrapper::ElasticOrthotropicThreeDimensionalWrapper(int tag, double Ex, double Ey,
	double Ez, double vxy, double vyz, double vzx, double Gxy, double Gyz,
	double Gzx, double rho) {
	_NDMaterial = new ElasticOrthotropicThreeDimensional(tag, Ex, Ey, Ez, vxy, vyz, vzx, Gxy, Gyz, Gzx, rho);
}

FSAMWrapper::FSAMWrapper(int tag,
	double RHO,
	UniaxialMaterialWrapper^ s1,
	UniaxialMaterialWrapper^ s2,
	UniaxialMaterialWrapper^ c1,
	UniaxialMaterialWrapper^ c2,
	UniaxialMaterialWrapper^ cA1,
	UniaxialMaterialWrapper^ cA2,
	UniaxialMaterialWrapper^ cB1,
	UniaxialMaterialWrapper^ cB2,
	double ROUX,
	double ROUY,
	double NU,
	double ALFADOW) {
	_NDMaterial = new FSAM(tag, RHO, s1->_UniaxialMaterial, s2->_UniaxialMaterial, c1->_UniaxialMaterial, c2->_UniaxialMaterial, cA1->_UniaxialMaterial,
		cA2->_UniaxialMaterial, cB1->_UniaxialMaterial, cB2->_UniaxialMaterial, ROUX, ROUY, NU, ALFADOW);
}

InitialStateAnalysisWrapperWrapper::InitialStateAnalysisWrapperWrapper(int tag, NDMaterialWrapper^ mainMat, int ndim) {
	_NDMaterial = new InitialStateAnalysisWrapper(tag, *mainMat->_NDMaterial, ndim);
}

InitStressNDMaterialWrapper::InitStressNDMaterialWrapper(int tag, NDMaterialWrapper^ material, VectorWrapper^ sigInit, int ndim) {
	_NDMaterial = new InitStressNDMaterial(tag, *material->_NDMaterial, *sigInit->_Vector, ndim);
}

LinearCapWrapper::LinearCapWrapper(int    tag,
	double G,
	double K,
	double rho,
	double theta,
	double alpha,
	double T,
	int ndm,
	double pTol_k) {
	_NDMaterial = new LinearCap(tag, G, K, rho, theta, alpha, T, ndm, pTol_k);
}

ManzariDafaliasWrapper::ManzariDafaliasWrapper(int tag, int classTag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double z_max, double cz, double mDen,
	int integrationScheme, int tangentType, int JacoType, double TolF, double TolR) {
	_NDMaterial = new ManzariDafalias(tag, classTag, G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, mDen, integrationScheme, tangentType, JacoType, TolF, TolR);
}

ManzariDafalias3DWrapper::ManzariDafalias3DWrapper(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double z_max, double cz, double mDen, int integrationScheme,
	int tangentType, int JacoType, double TolF, double TolR) {
	_NDMaterial = new ManzariDafalias3D(tag, G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, mDen, integrationScheme, tangentType, JacoType, TolF, TolR);
}

ManzariDafalias3DROWrapper::ManzariDafalias3DROWrapper(int tag, double G0, double nu, double B, double a1, double gamma1, double e_init, double Mc, double c,
	double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd,
	double z_max, double cz, double mDen, double kappa, int integrationScheme, int tangentType, int JacoType,
	double TolF, double TolR) {
	_NDMaterial = new ManzariDafalias3DRO(tag, G0, nu, B, a1, gamma1, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, mDen, kappa, integrationScheme, tangentType, JacoType, TolF, TolR);
}

ManzariDafaliasPlaneStrainWrapper::ManzariDafaliasPlaneStrainWrapper(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double z_max, double cz, double mDen, int integrationScheme,
	int tangentType, int JacoType, double TolF, double TolR) {
	_NDMaterial = new ManzariDafaliasPlaneStrain(tag, G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, mDen, integrationScheme, tangentType, JacoType, TolF, TolR);
}

ManzariDafaliasPlaneStrainROWrapper::ManzariDafaliasPlaneStrainROWrapper(int tag, double G0, double nu, double B, double a1, double gamma1, double e_init, double Mc, double c,
	double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd,
	double z_max, double cz, double mDen, double kappa, int integrationScheme, int tangentType, int JacoType,
	double TolF, double TolR) {
	_NDMaterial = new ManzariDafaliasPlaneStrainRO(tag, G0, nu, B, a1, gamma1, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, mDen, kappa, integrationScheme, tangentType, JacoType, TolF, TolR);
}

ManzariDafaliasROWrapper::ManzariDafaliasROWrapper(int tag, int classTag, double G0, double nu, double B, double a1, double gamma1, double e_init, double Mc, double c,
	double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd,
	double z_max, double cz, double mDen, double kappa, int integrationScheme, int tangentType, int JacoType,
	double TolF, double TolR) {
	_NDMaterial = new ManzariDafaliasRO(tag, classTag, G0, nu, B, a1, gamma1, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, mDen, kappa, integrationScheme, tangentType, JacoType, TolF, TolR);
}

PlasticDamageConcrete3dWrapper::PlasticDamageConcrete3dWrapper(int tag,
	double E,
	double nu,
	double ft,
	double fc,
	double beta,
	double Ap,
	double An,
	double Bn) {
	_NDMaterial = new PlasticDamageConcrete3d(tag, E, nu, ft, fc, beta, Ap, An, Bn);
}

PlasticDamageConcretePlaneStressWrapper::PlasticDamageConcretePlaneStressWrapper(int tag,
	double E,
	double nu,
	double ft,
	double fc,
	double beta,
	double Ap,
	double An,
	double Bn) {
	_NDMaterial = new PlasticDamageConcretePlaneStress(tag, E, nu, ft, fc, beta, Ap, An, Bn);
}

PM4SandWrapper::PM4SandWrapper(int tag, int classTag, double Dr, double G0, double hp0, double mDen, double P_atm, double h0, double emax,
	double emin, double nb, double nd, double Ado, double z_max, double cz,
	double ce, double phi_cv, double nu, double Cgd, double Cdr, double Ckaf, double Q,
	double R, double m, double Fsed_min, double p_sdeo, int integrationScheme, int tangentType,
	double TolF, double TolR) {
	_NDMaterial = new PM4Sand(tag, classTag, Dr, G0, hp0, mDen, P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ce, phi_cv, nu, Cgd, Cdr, Ckaf, Q, R, m, Fsed_min, p_sdeo, integrationScheme, tangentType, TolF, TolR);
}

PM4SandWrapper::PM4SandWrapper(int tag, double Dr, double G0, double hp0, double mDen, double P_atm, double h0, double emax,
	double emin, double nb, double nd, double Ado, double z_max, double cz,
	double ce, double phi_cv, double nu, double Cgd, double Cdr, double Ckaf, double Q,
	double R, double m, double Fsed_min, double p_sdeo, int integrationScheme, int tangentType,
	double TolF, double TolR) {
	_NDMaterial = new PM4Sand(tag, Dr, G0, hp0, mDen, P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ce, phi_cv, nu, Cgd, Cdr, Ckaf, Q, R, m, Fsed_min, p_sdeo, integrationScheme, tangentType, TolF, TolR);
}

PM4SiltWrapper::PM4SiltWrapper(int tag, int classTag, double Su, double Su_rate, double G0, double hpo, double mDen, double Fsu , double P_atm , double nu ,
	double nG , double h0 , double einit , double lambda , double phi_cv , double nbwet , double nbdry , double nd ,
	double Ado , double ru_max , double z_max , double cz , double ce, double Cgd , double Ckaf , double m ,
	double CG_consol , int integrationScheme , int tangentType , double TolF , double TolR ) {
	_NDMaterial = new PM4Silt( tag, classTag,  Su,  Su_rate,  G0,  hpo,  mDen,  Fsu,  P_atm,  nu,
		 nG,  h0,  einit,  lambda,  phi_cv,  nbwet,  nbdry,  nd,
		 Ado,  ru_max,  z_max,  cz,  ce,  Cgd,  Ckaf,  m,
		 CG_consol,  integrationScheme,  tangentType,  TolF,  TolR);
}

PM4SiltWrapper::PM4SiltWrapper(int tag, double Su, double Su_rate, double G0, double hpo, double mDen, double Fsu, double P_atm, double nu,
	double nG, double h0, double einit, double lambda, double phi_cv, double nbwet, double nbdry, double nd,
	double Ado, double ru_max, double z_max, double cz, double ce, double Cgd, double Ckaf, double m,
	double CG_consol, int integrationScheme, int tangentType, double TolF, double TolR) {
	_NDMaterial = new PM4Silt(tag, Su, Su_rate, G0, hpo, mDen, Fsu, P_atm, nu,
		nG, h0, einit, lambda, phi_cv, nbwet, nbdry, nd,
		Ado, ru_max, z_max, cz, ce, Cgd, Ckaf, m,
		CG_consol, integrationScheme, tangentType, TolF, TolR);
}

SimplifiedJ2Wrapper::SimplifiedJ2Wrapper(int tag,
	int nd,
	double G,
	double K,
	double sigmaY0,
	double H_kin,
	double H_iso) {
	_NDMaterial = new SimplifiedJ2(tag, nd, G, K, sigmaY0, H_kin, H_iso);
}

DruckerPragerThermalWrapper::DruckerPragerThermalWrapper(int tag, int classTag, double bulk, double shear,
	double s_y, double r, double r_bar, double Kinfinity, double Kinit,
	double d1, double d2, double H, double t, double massDen, double atm) {
	_NDMaterial = new DruckerPragerThermal(tag, classTag, bulk, shear, s_y, r, r_bar, Kinfinity, Kinit, d1, d2, H, t, massDen, atm);
}

DruckerPrager3DWrapper::DruckerPrager3DWrapper(int tag, double bulk, double shear,
	double s_y, double r, double r_bar, double Kinfinity, double Kinit,
	double d1, double d2, double H, double t, double massDen, double atm) {
	_NDMaterial = new DruckerPrager3D(tag, bulk, shear, s_y, r, r_bar, Kinfinity, Kinit, d1, d2, H, t, massDen, atm);
}

BoundingCamClayPlaneStrainWrapper::BoundingCamClayPlaneStrainWrapper(int tag, double mDen, double c, double bulk, double OCR, double mu_o,
	double alpha, double lambda, double h, double m) {
	_NDMaterial = new BoundingCamClayPlaneStrain(tag, mDen, c, bulk, OCR, mu_o, alpha, lambda, h, m);
}

BoundingCamClay3DWrapper::BoundingCamClay3DWrapper(int tag, double mDen, double c, double bulk, double OCR, double mu_o,
	double alpha, double lambda, double h, double m) {
	_NDMaterial = new BoundingCamClay3D(tag, mDen, c, bulk, OCR, mu_o, alpha, lambda, h, m);
}

BoundingCamClayWrapper::BoundingCamClayWrapper(int tag, int classTag, double massDen, double C, double bulk, double OCR,
	double mu_o, double Alpha, double lambda, double h, double m) {
	_NDMaterial = new BoundingCamClay(tag, classTag, massDen, C, bulk, OCR, mu_o, Alpha, lambda, h, m);
}

AcousticMediumWrapper::AcousticMediumWrapper(int tag, int classTag, double k, double rho, double gamma) {
	_NDMaterial = new AcousticMedium(tag, classTag, k, rho, gamma);
}

ExternalNDMaterialWrapper::ExternalNDMaterialWrapper(int tag)
{
	_NDMaterial = new ExternalNDMaterial(tag);
}
