/***************************************************************************
                          CombinedIsoKin2D02.cpp  -  description
                             -------------------
    begin                : Fri Jul 12 2002
    email                : rkaul@ce-blume215-pent-2.stanford.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *                                                                         *
 ***************************************************************************/

#include "CombinedIsoKin2D02.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CombinedIsoKin2D02::CombinedIsoKin2D02(int tag,  double min_iso_factor,
				double iso_ratio, double kin_ratio,
				YieldSurface_BC  &lim_surface,
				PlasticHardeningMaterial &kinX,
				PlasticHardeningMaterial &kinY,
				PlasticHardeningMaterial &isoXPos,
				PlasticHardeningMaterial &isoXNeg,
				PlasticHardeningMaterial &isoYPos,
				PlasticHardeningMaterial &isoYNeg,
				bool isDeformable,
				int  algo, double resfact, double appfact, double dir)
:BkStressLimSurface2D(tag, -1, min_iso_factor,
				iso_ratio, kin_ratio,
				lim_surface, kinX, kinY,
				isoXPos, isoXNeg, isoYPos, isoYNeg, algo, resfact, appfact, dir)
{
	deformable = isDeformable;
}

CombinedIsoKin2D02::~CombinedIsoKin2D02()
{
}


YS_Evolution * CombinedIsoKin2D02::getCopy()
{
	CombinedIsoKin2D02 *theCopy = new 
					CombinedIsoKin2D02(
					  getTag(), minIsoFactor,
					  isotropicRatio_orig, kinematicRatio_orig,
					  *limSurface, *kinMatX, *kinMatY,
					  *isoMatXPos, *isoMatXNeg, *isoMatYPos, *isoMatYNeg,
					  deformable, resAlgo, resFactor, appFactor, direction_orig);
	return theCopy;
}


void CombinedIsoKin2D02::Print(OPS_Stream & s, int flag)
{
	s << "CombinedIsoKin2D02 \n";
	s << "iso_Ratio = " << isotropicRatio << "\n";
	s << "isotropicFactor_hist = " << isotropicFactor_hist;
	s << "translateX       = " << translate(0) << ",\ttranslateY = " << translate(1) << "\n";
	s << "\n";

}
