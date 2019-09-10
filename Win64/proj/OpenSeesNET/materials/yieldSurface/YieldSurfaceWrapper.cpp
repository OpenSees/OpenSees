#include "stdafx.h"
#include "YieldSurfaceWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Materials::YieldSurfaces;


ExponReducingWrapper::ExponReducingWrapper(int tag, double kp0, double alfa)
{
	_PlasticHardeningMaterial = new ExponReducing(tag, kp0, alfa);
}

ExponReducingWrapper::ExponReducingWrapper(int tag, double kp0, double alfa, double res_fact)
{
	_PlasticHardeningMaterial = new ExponReducing(tag, kp0, alfa, res_fact);
}

MultiLinearKpWrapper::MultiLinearKpWrapper(int tag, VectorWrapper^ sum_plas_defo, VectorWrapper^ kp)
{
	_PlasticHardeningMaterial = new MultiLinearKp(tag, *sum_plas_defo->_Vector, *kp->_Vector);
}

NullPlasticMaterialWrapper::NullPlasticMaterialWrapper(int tag)
{
	_PlasticHardeningMaterial = new NullPlasticMaterial(tag);
}

Attalla2DWrapper::Attalla2DWrapper(int tag, double xmax, double ymax,
	YS_EvolutionWrapper^ model)
{
	_YieldSurface_BC2D = new Attalla2D(tag, xmax, ymax,
		*model->_YS_Evolution);
}

Attalla2DWrapper::Attalla2DWrapper(int tag, double xmax, double ymax,
	YS_EvolutionWrapper^ model, double a01, double a02, double a03,
	double a04, double a05, double a06)
{
	_YieldSurface_BC2D = new Attalla2D(tag, xmax, ymax,
		*model->_YS_Evolution);
}

ElTawil2DWrapper::ElTawil2DWrapper(int tag, double xbal, double ybal, double ypos, double yneg,
	YS_EvolutionWrapper^ model)
{
	_YieldSurface_BC2D = new ElTawil2D( tag,  xbal,  ybal,  ypos,  yneg,
		*model->_YS_Evolution);
}

ElTawil2DWrapper::ElTawil2DWrapper(int tag, double xbal, double ybal, double ypos, double yneg,
	YS_EvolutionWrapper^ model, double cz, double ty)
{
	_YieldSurface_BC2D = new ElTawil2D(tag, xbal, ybal, ypos, yneg,
		*model->_YS_Evolution, cz, ty);
}

ElTawil2DUnSymWrapper::ElTawil2DUnSymWrapper(int tag, double xPosBal, double yPosBal,
	double xNegBal, double yNegBal,
	double ypos, double yneg,
	YS_EvolutionWrapper^ model)
{
	_YieldSurface_BC2D = new ElTawil2DUnSym( tag,  xPosBal,  yPosBal,
		 xNegBal,  yNegBal,
		 ypos,  yneg,
		*model->_YS_Evolution);
}

ElTawil2DUnSymWrapper::ElTawil2DUnSymWrapper(int tag, double xPosBal, double yPosBal,
	double xNegBal, double yNegBal,
	double ypos, double yneg,
	YS_EvolutionWrapper^ model, double czPos, double tyPos,
	double czNeg, double tyNeg)
{
	_YieldSurface_BC2D = new ElTawil2DUnSym(tag, xPosBal, yPosBal,
		xNegBal, yNegBal,
		ypos, yneg,
		*model->_YS_Evolution, czPos, tyPos,
		czNeg, tyNeg);
}


Hajjar2DWrapper::Hajjar2DWrapper(int tag, double xmax, double ymax,
	YS_EvolutionWrapper^ model,
	double centroid_y, double c1, double c2, double c3) {
	_YieldSurface_BC2D = new Hajjar2D(tag, xmax, ymax,
		*model->_YS_Evolution,
		centroid_y, c1, c2, c3);
}
Hajjar2DWrapper::Hajjar2DWrapper(int tag, YS_EvolutionWrapper^ model,
	double D, double b, double t, double fc_, double fy_) {
	_YieldSurface_BC2D = new Hajjar2D(tag,
		*model->_YS_Evolution,
		D, b, t, fc_, fy_);
}

NullYS2DWrapper::NullYS2DWrapper(int tag) {
	_YieldSurface_BC2D = new NullYS2D(tag);
}

Orbison2DWrapper::Orbison2DWrapper(int tag, double xmax, double ymax,  YS_EvolutionWrapper^ model) {
	_YieldSurface_BC2D = new Orbison2D(tag, xmax, ymax,
		*model->_YS_Evolution);
}

CombinedIsoKin2D01Wrapper::CombinedIsoKin2D01Wrapper(int tag,
	double iso_ratio, double kin_ratio,
	double shr_iso_ratio, double shr_kin_ratio,
	double min_iso_factor,
	PlasticHardeningMaterialWrapper^ kpx_pos,
	PlasticHardeningMaterialWrapper^  kpx_neg,
	PlasticHardeningMaterialWrapper^ kpy_pos,
	PlasticHardeningMaterialWrapper^  kpy_neg,
	bool isDeformable, double dir) {
	_YS_Evolution2D = new CombinedIsoKin2D01(tag,
		iso_ratio, kin_ratio,
		shr_iso_ratio, shr_kin_ratio,
		min_iso_factor,
		*kpx_pos->_PlasticHardeningMaterial,
		*kpx_neg->_PlasticHardeningMaterial,
		*kpy_pos->_PlasticHardeningMaterial,
		*kpy_neg->_PlasticHardeningMaterial,
		isDeformable, dir);
}

CombinedIsoKin2D02Wrapper::CombinedIsoKin2D02Wrapper(int tag, double min_iso_factor,
	double iso_ratio, double kin_ratio,
	YieldSurface_BCWrapper^  lim_surface,
	PlasticHardeningMaterialWrapper^ kinX,
	PlasticHardeningMaterialWrapper^  kinY,
	PlasticHardeningMaterialWrapper^ isoXPos,
	PlasticHardeningMaterialWrapper^  isoXNeg,
	PlasticHardeningMaterialWrapper^ isoYPos,
	PlasticHardeningMaterialWrapper^  isoYNeg,
	bool isDeformable,
	int  algo, double resfact, double appfact, double dir) {
	_YS_Evolution2D = new CombinedIsoKin2D02(tag, min_iso_factor,
		iso_ratio, kin_ratio,
		*lim_surface->_YieldSurface_BC,
		*kinX->_PlasticHardeningMaterial,
		*kinY->_PlasticHardeningMaterial,
		*isoXPos->_PlasticHardeningMaterial,
		*isoXNeg->_PlasticHardeningMaterial,
		*isoYPos->_PlasticHardeningMaterial,
		*isoYNeg->_PlasticHardeningMaterial,
		isDeformable,
		algo, resfact, appfact, dir);
}

Isotropic2D01Wrapper::Isotropic2D01Wrapper(int tag, double min_iso_factor,
	PlasticHardeningMaterialWrapper^ kpx,
	PlasticHardeningMaterialWrapper^  kpy) {

	_YS_Evolution2D = new Isotropic2D01(tag, min_iso_factor,
		*kpx->_PlasticHardeningMaterial,
		*kpy->_PlasticHardeningMaterial);
}

Kinematic2D01Wrapper::Kinematic2D01Wrapper(int tag, double min_iso_factor,
	PlasticHardeningMaterialWrapper^ kpx,
	PlasticHardeningMaterialWrapper^  kpy, double dir) {

	_YS_Evolution2D = new Kinematic2D01(tag, min_iso_factor,
		*kpx->_PlasticHardeningMaterial,
		*kpy->_PlasticHardeningMaterial, dir);
}

Kinematic2D02Wrapper::Kinematic2D02Wrapper(int tag, double min_iso_factor,
	YieldSurface_BCWrapper^ lim_surface,
	PlasticHardeningMaterialWrapper^ kpx,
	PlasticHardeningMaterialWrapper^  kpy, 
	int algo, double resfact, double appfact, double dir) {

	_YS_Evolution2D = new Kinematic2D02(tag, min_iso_factor,
		*lim_surface->_YieldSurface_BC,
		*kpx->_PlasticHardeningMaterial,
		*kpy->_PlasticHardeningMaterial, algo, resfact, appfact, dir);
}

NullEvolutionWrapper::NullEvolutionWrapper(int tag, double isox) {

	_YS_Evolution = new NullEvolution(tag, isox);
}

NullEvolutionWrapper::NullEvolutionWrapper(int tag, double isox, double isoy) {

	_YS_Evolution = new NullEvolution(tag, isox, isoy);
}

NullEvolutionWrapper::NullEvolutionWrapper(int tag, double isox, double isoy, double isoz) {

	_YS_Evolution = new NullEvolution(tag, isox, isoy, isoz);
}

PeakOriented2D01Wrapper::PeakOriented2D01Wrapper(int tag, double min_iso_factor,
	PlasticHardeningMaterialWrapper^ kpx,
	PlasticHardeningMaterialWrapper^  kpy) {

	_YS_Evolution2D = new PeakOriented2D01(tag, min_iso_factor,
		*kpx->_PlasticHardeningMaterial,
		*kpy->_PlasticHardeningMaterial);
}

PeakOriented2D02Wrapper::PeakOriented2D02Wrapper(int tag, double min_iso_factor,
	YieldSurface_BCWrapper^ lim_surface,
	PlasticHardeningMaterialWrapper^ kinX,
	PlasticHardeningMaterialWrapper^  kinY,
	PlasticHardeningMaterialWrapper^ isoX,
	PlasticHardeningMaterialWrapper^  isoY, int algo) {

	_YS_Evolution2D = new PeakOriented2D02(tag, min_iso_factor,
		*lim_surface->_YieldSurface_BC,
		*kinX->_PlasticHardeningMaterial,
		*kinY->_PlasticHardeningMaterial,
		*isoX->_PlasticHardeningMaterial,
		*isoY->_PlasticHardeningMaterial,
		algo);
}