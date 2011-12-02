/***************************************************************************
                          PeakOriented2D02.cpp  -  description
                             -------------------
    begin                : Fri Jul 12 2002
    email                : rkaul@ce-blume215-pent-2.stanford.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *                                                                         *
 ***************************************************************************/

#include "PeakOriented2D02.h"
NullPlasticMaterial PeakOriented2D02::nullMat(-1);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PeakOriented2D02::PeakOriented2D02(int tag, double min_iso_factor,
					YieldSurface_BC &lim_surface,
					PlasticHardeningMaterial &kinX,
					PlasticHardeningMaterial &kinY,
					PlasticHardeningMaterial &isoX,
					PlasticHardeningMaterial &isoY,
					int algo)
:BkStressLimSurface2D(tag, -1, min_iso_factor, 0.5, 0.5,
                      lim_surface,
                      kinX, kinY,
                      isoX, isoX, isoY, isoY, algo, 0.0, 0.0, 0.0)
{
                                                  
}

PeakOriented2D02::~PeakOriented2D02()
{
	
}

YS_Evolution * PeakOriented2D02::getCopy()
{
 PeakOriented2D02 *theMat = new
           PeakOriented2D02( getTag(), minIsoFactor, *limSurface, *kinMatX, *kinMatY,
                              *isoMatXPos, *isoMatYPos, resAlgo);
                              
	return theMat;
}

void PeakOriented2D02::Print(OPS_Stream & s, int flag)
{
	s << "PeakOriented2D02 \n";
	s << "iso_Ratio = " << isotropicRatio << "\n";
	s << "isotropicFactor_hist = " << isotropicFactor_hist;
	s << "translateX       = " << translate(0) << ",\ttranslateY = " << translate(1) << "\n";
	s << "\n";

}
